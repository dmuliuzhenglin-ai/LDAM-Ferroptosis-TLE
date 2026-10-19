# ============================================================================
# Script      : 05_immune_GSVA_enrichment.R
# Purpose     : GO/KEGG/GSEA enrichment, ssGSEA immune infiltration and sample-level GSVA LDAM/ferroptosis scoring in bulk data.
# Manuscript  : Methods 2.2 / Results 3.3
# Input       : Candidate gene lists from 04; expression matrices; genesets/
# Output      : results/05_* enrichment and correlation tables/figures
# Dependencies: GSVA, clusterProfiler, data.table, dplyr, enrichplot, ggplot2, ggpubr, org.Hs.eg.db, pheatmap
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

# ==============================================================================
# 05_immune_GSVA_enrichment.R
# ① 候选基因 GO/KEGG/GSEA 富集（clusterProfiler）
# ② ssGSEA 免疫浸润（含小胶质细胞/巨噬细胞）+ hub 基因相关性
# ③ GSVA 计算每样本 LDAM/铁死亡评分并检验相关性（bulk 层面，与单细胞互证）
# 依赖：先运行 00 / 01 / 02 / 04
# ==============================================================================

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results"); FIG_DIR <- file.path(BASE_DIR, "figures")
GS_DIR   <- file.path(BASE_DIR, "genesets")

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(GSVA); library(ggplot2); library(ggpubr); library(pheatmap)
})

LDAM_genes <- readLines(file.path(GS_DIR, "LDAM_geneset.txt"))
Ferroptosis_genes <- readLines(file.path(GS_DIR, "Ferroptosis_geneset.txt"))
hub_genes <- readLines(file.path(OUT_DIR, "04_hub_genes.txt"))

## ---------- 5.1 GO / KEGG 富集（候选基因） ----------
candidates <- readLines(file.path(OUT_DIR, "04_candidate_genes.txt"))
entrez <- bitr(candidates, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ekegg <- enrichKEGG(entrez$ENTREZID, organism = "hsa", pAdjustMethod = "BH",
                    qvalueCutoff = 0.05)
fwrite(as.data.frame(ego), file.path(OUT_DIR, "05_GO_BP_candidates.csv"))
fwrite(as.data.frame(ekegg), file.path(OUT_DIR, "05_KEGG_candidates.csv"))

pdf(file.path(FIG_DIR, "05_GO_BP_dotplot.pdf"), width = 9, height = 7)
print(dotplot(ego, showCategory = 15)); dev.off()
if (nrow(as.data.frame(ekegg)) > 0) {
  pdf(file.path(FIG_DIR, "05_KEGG_dotplot.pdf"), width = 9, height = 7)
  print(dotplot(ekegg, showCategory = 15)); dev.off()
}

## ---------- 5.2 GSEA（GSE168375 全基因排序） ----------
de2 <- fread(file.path(OUT_DIR, "02_GSE168375_DEG_Epilepsy_vs_Control.csv"))
gene_list <- de2$log2FoldChange; names(gene_list) <- de2$gene
gene_list <- sort(gene_list, decreasing = TRUE); gene_list <- gene_list[!duplicated(names(gene_list))]
entrez_list <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
gl <- gene_list[entrez_list$SYMBOL]; names(gl) <- entrez_list$ENTREZID
gsea_res <- gseKEGG(gl, organism = "hsa", pvalueCutoff = 0.05)
fwrite(as.data.frame(gsea_res), file.path(OUT_DIR, "05_GSEA_KEGG_GSE168375.csv"))
if (nrow(as.data.frame(gsea_res)) > 0) {
  pdf(file.path(FIG_DIR, "05_GSEA_KEGG_ridge.pdf"), width = 10, height = 7)
  print(ridgeplot(gsea_res, showCategory = 15)); dev.off()
}

## ---------- 5.3 ssGSEA 免疫浸润（GSE168375） ----------
# 28 种免疫细胞 marker（Charoentong 2017 简化版，含小胶质细胞/巨噬细胞）
immune_sigs <- list(
  Microglia   = c("P2RY12","TMEM119","CX3CR1","CSF1R","AIF1","C1QA","C1QB","C1QC","GPR34","SIGLEC8"),
  Macrophage  = c("CD68","CD163","MRC1","MSR1","MARCO","FCGR3A","LYZ","C1QA"),
  Monocyte    = c("CD14","FCGR3A","S100A8","S100A9","VCAN","LYZ"),
  T_cell_CD4  = c("CD3D","CD3E","CD4","IL7R","LTB"),
  T_cell_CD8  = c("CD3D","CD3E","CD8A","CD8B","GZMB","NKG7"),
  B_cell      = c("CD19","MS4A1","CD79A","CD79B"),
  NK_cell     = c("NCAM1","KLRD1","KLRB1","GNLY","NKG7"),
  Neutrophil  = c("FCGR3B","CSF3R","S100A8","S100A9","CXCR2"),
  DC          = c("ITGAX","HLA-DRA","HLA-DRB1","CLEC10A","CD1C"),
  Mast_cell   = c("TPSAB1","TPSB2","CPA3","KIT")
)

expr_vst <- fread(file.path(OUT_DIR, "01_GSE134697_vst_normalized.csv")) |>
  tibble::column_to_rownames("gene") |> as.matrix()
meta <- fread(file.path(OUT_DIR, "01_GSE134697_sample_meta.csv"))
meta$group <- meta$tissue   # Hippocampus vs Neocortex

immune_score <- gsva(gsvaParam(expr_vst, immune_sigs))   # 细胞类型 x 样本
fwrite(data.frame(cell_type = rownames(immune_score), immune_score),
       file.path(OUT_DIR, "05_ssGSEA_immune_scores.csv"))

# 组间差异箱线图（重点：Microglia / Macrophage）
imm_long <- as.data.frame(t(immune_score)) |>
  tibble::rownames_to_column("sample") |>
  left_join(meta, by = "sample") |>
  tidyr::pivot_longer(cols = -any_of(c("sample", "tissue", "patient", "group")),
                      names_to = "cell", values_to = "score")
p_imm <- ggplot(imm_long, aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~cell, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("Neocortex" = "#7BA7BC", "Hippocampus" = "#C0392B")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3) +
  labs(title = "GSE134697 ssGSEA 免疫浸润：TLE 海马 vs 自体皮层") + theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(FIG_DIR, "05_immune_infiltration_boxplot.pdf"), p_imm,
       width = 12, height = 6)

## ---------- 5.4 hub 基因表达 vs 小胶质细胞/巨噬细胞浸润相关性 ----------
mg_score <- immune_score["Microglia", ]
mp_score <- immune_score["Macrophage", ]
cor_res <- lapply(hub_genes, function(g) {
  if (!g %in% rownames(expr_vst)) return(NULL)
  data.frame(gene = g,
             Microglia_cor = cor(expr_vst[g, ], mg_score, method = "spearman"),
             Microglia_p   = cor.test(expr_vst[g, ], mg_score, method = "spearman")$p.value,
             Macrophage_cor = cor(expr_vst[g, ], mp_score, method = "spearman"),
             Macrophage_p   = cor.test(expr_vst[g, ], mp_score, method = "spearman")$p.value)
}) |> bind_rows()
print(cor_res)
fwrite(cor_res, file.path(OUT_DIR, "05_hub_vs_microglia_macrophage_cor.csv"))

## ---------- 5.5 GSVA：LDAM 评分 vs 铁死亡评分相关性（bulk 层面互证） ----------
gs2 <- list(LDAM = LDAM_genes, Ferroptosis = Ferroptosis_genes)
score2 <- gsva(gsvaParam(expr_vst, gs2))

cor_test <- cor.test(score2["LDAM", ], score2["Ferroptosis", ], method = "spearman")
df_s <- data.frame(LDAM = score2["LDAM", ], Ferroptosis = score2["Ferroptosis", ],
                   group = meta$group[match(colnames(score2), meta$sample)])
p_cor <- ggplot(df_s, aes(LDAM, Ferroptosis, color = group)) +
  geom_point(size = 2.5) + geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  scale_color_manual(values = c("Neocortex" = "#7BA7BC", "Hippocampus" = "#C0392B")) +
  annotate("text", x = min(df_s$LDAM), y = max(df_s$Ferroptosis), hjust = 0,
           label = paste0("Spearman r = ", round(cor_test$estimate, 3),
                          ", P = ", signif(cor_test$p.value, 3))) +
  labs(title = "GSE134697：LDAM 评分 vs 铁死亡评分（GSVA）") + theme_bw()
ggsave(file.path(FIG_DIR, "05_LDam_Ferroptosis_correlation_bulk.pdf"),
       p_cor, width = 7, height = 6)

# 组间评分差异
p_score <- df_s |>
  tidyr::pivot_longer(c(LDAM, Ferroptosis), names_to = "pathway", values_to = "score") |>
  ggplot(aes(group, score, fill = group)) +
  geom_boxplot() + facet_wrap(~pathway, scales = "free") +
  scale_fill_manual(values = c("Neocortex" = "#7BA7BC", "Hippocampus" = "#C0392B")) +
  stat_compare_means(method = "wilcox.test") + theme_bw()
ggsave(file.path(FIG_DIR, "05_LDam_Ferroptosis_score_boxplot.pdf"),
       p_score, width = 8, height = 5)

message("bulk 层面 LDAM-铁死亡相关性: r = ", round(cor_test$estimate, 3),
        ", P = ", signif(cor_test$p.value, 3))
message("=== 05_immune_GSVA_enrichment.R 运行完毕 ===")
