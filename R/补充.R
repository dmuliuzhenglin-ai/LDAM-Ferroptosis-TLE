# ============================================================================
# Script      : 补充.R
# Purpose     : Supplementary single-cell analyses: LDAM-high vs -low DE + GO, hub-gene cell-type specificity, microglia pseudobulk DE, hub expression in mouse KA.
# Manuscript  : Methods 2.4 supplementary
# Input       : results/06 rds; mouse KA data; genesets/
# Output      : results/09_* supplementary tables/figures
# Dependencies: DESeq2, Seurat, clusterProfiler, data.table, dplyr, ggplot2, org.Hs.eg.db, pROC
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 09_补充分析.R（最终版：所有已知问题已修复）
## 补1：LDAM-high vs LDAM-low 差异基因 + GO
## 补2：hub 基因细胞类型特异性
## 补3：小胶质细胞 pseudobulk 差异分析
## 补4：hub 基因在小鼠 KA 模型中的表达
## ============================================================

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results"); FIG_DIR <- file.path(BASE_DIR, "figures")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(data.table)
})

## ---------- 补4：hub 基因在小鼠 KA 模型中的表达 ----------
de_m <- fread(file.path(OUT_DIR, "08_GSE88992_mouse_DEG_KA_vs_Saline.csv"))
hub_mouse <- de_m[tolower(gene) %in% tolower(c("LPCAT2","S100A10","CD9","SYNGR2"))]
cat("== 补4：hub 基因在小鼠 KA 模型海马中的表达 ==\n")
print(hub_mouse[, .(gene, logFC, P.Value, adj.P.Val)])

## ---------- 读入 06 的对象（立即合并层，防 Seurat v5 报错） ----------
mg  <- readRDS(file.path(OUT_DIR, "06_GSE190452_microglia.rds"))
obj <- readRDS(file.path(OUT_DIR, "06_GSE190452_annotated_full.rds"))
mg  <- JoinLayers(mg)
obj <- JoinLayers(obj)
hub_genes <- readLines(file.path(OUT_DIR, "04_hub_genes.txt"))

## ---------- 补1：LDAM-high vs LDAM-low 差异基因 + GO ----------
Idents(mg) <- factor(mg$LDAM_level, levels = c("LDAM-low", "LDAM-high"))
de_ldam <- FindMarkers(mg, ident.1 = "LDAM-high", ident.2 = "LDAM-low",
                       min.pct = 0.1, logfc.threshold = 0.25)
de_ldam$gene <- rownames(de_ldam)
fwrite(de_ldam, file.path(OUT_DIR, "09_DEG_LDAMhigh_vs_LDAMlow.csv"))
cat("== 补1：LDAM-high vs low 差异基因数（padj<0.05）==\n")
print(table(de_ldam$p_val_adj < 0.05))
cat("LDAM-high 上调前 20：\n")
# [debug] print(head(de_ldam[order(-de_ldam$avg_log2FC), ], 20)[, c("gene","avg_log2FC","p_val_adj")])

sig_up <- de_ldam$gene[de_ldam$p_val_adj < 0.05 & de_ldam$avg_log2FC > 0]
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE) && length(sig_up) >= 10) {
  library(clusterProfiler); library(org.Hs.eg.db)
  ego <- enrichGO(sig_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pAdjustMethod = "BH", qvalueCutoff = 0.05)
  fwrite(as.data.frame(ego), file.path(OUT_DIR, "09_GO_LDAMhigh_up_genes.csv"))
  if (nrow(as.data.frame(ego)) > 0) {
    p_go <- dotplot(ego, showCategory = 15) +
      ggtitle("LDAM-high microglia upregulated genes - GO BP")
    ggsave(file.path(FIG_DIR, "09_GO_LDAMhigh.pdf"), p_go, width = 9, height = 7)
  }
  cat("GO 富集完成，显著通路数：", nrow(as.data.frame(ego)), "\n")
} else {
  message("clusterProfiler 未安装或上调基因太少，跳过 GO")
}

## ---------- 补2：hub 基因细胞类型特异性 ----------
hg <- intersect(hub_genes, rownames(obj))
if (length(hg) > 0) {
  p_hub_ct <- DotPlot(obj, features = hg, group.by = "cell_type") + RotatedAxis() +
    ggtitle("Hub genes expression across cell types")
  ggsave(file.path(FIG_DIR, "09_hub_celltype_specificity.pdf"),
         p_hub_ct, width = 6, height = 5)
  avg <- AverageExpression(obj, features = hg, group.by = "cell_type")$RNA
  cat("== 补2：hub 基因各细胞类型平均表达 ==\n")
  print(round(avg, 3))
}

## ---------- 补3：小胶质细胞 pseudobulk 差异分析 ----------
pb_counts <- AggregateExpression(mg, group.by = "orig.ident", assays = "RNA")$RNA
colnames(pb_counts) <- gsub("-", "_", colnames(pb_counts))   # 修复下划线被转成连字符
smeta <- unique(mg@meta.data[, c("orig.ident", "group")])
grp <- smeta$group[match(colnames(pb_counts), smeta$orig.ident)]
stopifnot(!any(is.na(grp)))
print(grp)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(round(as.matrix(pb_counts)),
                              data.frame(row.names = colnames(pb_counts),
                                         group = factor(grp, levels = c("Control", "TLE"))),
                              ~ group)
dds <- DESeq(dds)
res_mg <- as.data.frame(results(dds, contrast = c("group", "TLE", "Control")))
res_mg$gene <- rownames(res_mg)
fwrite(res_mg, file.path(OUT_DIR, "09_microglia_pseudobulk_DEG_TLE_vs_Control.csv"))
cat("== 补3：pseudobulk DEG 数（padj<0.05 且 |log2FC|>0.585）==\n")
print(table(res_mg$padj < 0.05 & abs(res_mg$log2FoldChange) > 0.585))
cat("hub 基因在 pseudobulk 中的表现：\n")
print(res_mg[res_mg$gene %in% hub_genes, c("gene","log2FoldChange","pvalue","padj")])

## hub 基因 pseudobulk 热图
vst_mg <- assay(vst(dds))
hg2 <- intersect(hub_genes, rownames(vst_mg))
if (length(hg2) > 0 && requireNamespace("pheatmap", quietly = TRUE)) {
  pdf(file.path(FIG_DIR, "09_hub_pseudobulk_heatmap.pdf"), width = 6, height = 4)
  pheatmap::pheatmap(vst_mg[hg2, ], scale = "row",
                     annotation_col = data.frame(group = grp, row.names = colnames(vst_mg)),
                     cluster_cols = FALSE, main = "Hub genes in microglia (pseudobulk)")
  dev.off()
}

cat("\n== 09 补充分析全部完成 ==\n")

library(data.table); library(pROC)
if (!exists("BASE")) BASE <- BASE_DIR  # alias of project root
OUT_DIR <- file.path(BASE, "results"); ADV_DIR <- file.path(BASE, "figures", "advanced")
hub_genes <- readLines(file.path(OUT_DIR, "04_hub_genes.txt"))
ex <- fread(file.path(OUT_DIR, "02_GSE168375_vst_normalized.csv"))
gcol <- names(ex)[1]; ex <- ex[!is.na(ex[[gcol]]) & ex[[gcol]] != ""]
exm <- as.matrix(ex[, -1, with = FALSE]); rownames(exm) <- ex[[gcol]]
mt <- fread(file.path(OUT_DIR, "02_GSE168375_sample_meta.csv"))
sm <- intersect(mt$sample, colnames(exm))
grp <- factor(ifelse(mt$group[match(sm, mt$sample)] %in% c("TLE","Epilepsy"), "Epilepsy", "Control"),
              levels = c("Control", "Epilepsy"))
df <- data.frame(grp = grp, t(exm[hub_genes, sm]))   # 43 行 x 5 列
str(df)

## LOOCV
pred <- numeric(nrow(df))
for (i in seq_len(nrow(df))) {
  fit <- glm(grp ~ ., data = df[-i, ], family = binomial)
  pred[i] <- predict(fit, newdata = df[i, , drop = FALSE], type = "response")
}
roc_loocv <- roc(df$grp, pred, quiet = TRUE); ci_l <- ci.auc(roc_loocv)
cat("LOOCV AUC =", round(auc(roc_loocv), 3), "(95%CI", round(ci_l[1],2), "-", round(ci_l[3],2), ")\n")

## 逐基因 ROC
roc_list <- lapply(hub_genes, function(g) roc(df$grp, df[[g]], quiet = TRUE)); names(roc_list) <- hub_genes
cols_g <- c("#C0392B","#5B8E7D","#D4A24E","#2C5F8A")
pdf(file.path(ADV_DIR, "11_3_LOOCV_pergene_ROC.pdf"), width = 11, height = 5)
par(mfrow = c(1, 2))
plot(roc_loocv, col = "#2C5F8A", lwd = 2, main = paste0("LOOCV ROC, AUC = ", round(auc(roc_loocv), 2)))
abline(0, 1, lty = 2, col = "grey60")
plot(roc_list[[1]], col = cols_g[1], lwd = 2, main = "Per-gene ROC")
for (k in 2:4) plot(roc_list[[k]], col = cols_g[k], lwd = 2, add = TRUE)
legend("bottomright", paste(hub_genes, round(sapply(roc_list, auc), 2)), col = cols_g, lwd = 2, cex = 0.8)
abline(0, 1, lty = 2, col = "grey60"); dev.off()

## Bootstrap CI（1000 次）
set.seed(42)
aucs <- replicate(1000, {
  idx <- sample(nrow(df), replace = TRUE)
  if (length(unique(df$grp[idx])) < 2) return(NA)
  fit <- tryCatch(glm(grp ~ ., data = df[idx, ], family = binomial), error = function(e) NULL)
  if (is.null(fit)) return(NA)
  auc(roc(df$grp[idx], predict(fit, type = "response"), quiet = TRUE))
})
aucs <- aucs[!is.na(aucs)]
cat("Bootstrap AUC =", round(mean(aucs), 3), "(95%CI",
    round(quantile(aucs, .025), 2), "-", round(quantile(aucs, .975), 2), ")\n")
fwrite(data.frame(metric = c("LOOCV_AUC","CI_lo","CI_hi","Boot_AUC","Boot_lo","Boot_hi"),
                  value = c(auc(roc_loocv), ci_l[1], ci_l[3], mean(aucs),
                            quantile(aucs,.025), quantile(aucs,.975))),
       file.path(OUT_DIR, "11_3_ROC_robustness.csv"))
cat("补3 完成\n")
