# ============================================================================
# Script      : 10.R
# Purpose     : Publication-grade main figures: assembles 7 multi-panel main figures (PDF + PNG) directly from results/.
# Manuscript  : Main Figures 1-7
# Input       : results/*.csv, *.rds from upstream scripts
# Output      : figures/main_figures/ (PDF + PNG)
# Dependencies: CellChat, Seurat, data.table, dplyr, ggplot2, pROC, patchwork, tidyr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 10_出版级图片.R
## 从 results/ 的 CSV 与 rds 直接出 7 张主图（PDF + PNG）
## 统一主题：Arial 风格、对照 #7BA7BC / 疾病 #C0392B
## 缺哪个文件就跳过对应 panel，不会中断
## ============================================================

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results"); FIG_DIR <- file.path(BASE_DIR, "figures")
MAIN_DIR <- file.path(FIG_DIR, "main_figures"); dir.create(MAIN_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(patchwork)
  library(pROC); library(tidyr)
})

## ---------- 统一主题与配色 ----------
theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
          plot.tag = element_text(face = "bold", size = base_size + 4),
          legend.position = "right",
          strip.background = element_rect(fill = "grey92", color = NA),
          strip.text = element_text(face = "bold"))
}
COL <- c(Control = "#7BA7BC", TLE = "#C0392B", Epilepsy = "#C0392B",
         Saline = "#7BA7BC", KA = "#C0392B",
         `LDAM-low` = "#7BA7BC", `LDAM-high` = "#C0392B")

save_fig <- function(p, name, w, h) {
  ggsave(file.path(MAIN_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(MAIN_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat(name, "已保存\n")
}

volcano_plot <- function(dt, title) {
  dt[, dir := fifelse(padj < 0.05 & log2FoldChange >  0.585, "Up",
                      fifelse(padj < 0.05 & log2FoldChange < -0.585, "Down", "NS"))]
  ggplot(dt, aes(log2FoldChange, -log10(pvalue), color = dir)) +
    geom_point(alpha = 0.55, size = 0.7) +
    scale_color_manual(values = c(Up = "#C0392B", Down = "#2C5F8A", NS = "grey85")) +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = 2, color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey60", linewidth = 0.3) +
    labs(title = title, x = "log2 fold change", y = "-log10 P") +
    theme_pub() + theme(legend.position = "none")
}

roc_from_csv <- function(expr_file, meta_file, genes, title) {
  ex <- fread(expr_file) |> tibble::column_to_rownames("gene") |> as.matrix()
  mt <- fread(meta_file)
  g <- intersect(genes, rownames(ex))
  sm <- intersect(mt$sample, colnames(ex))
  Xv <- t(ex[g, sm, drop = FALSE])
  grp <- factor(ifelse(mt$group[match(sm, mt$sample)] %in% c("TLE","Epilepsy","KA"), "Epilepsy", "Control"),
                levels = c("Control","Epilepsy"))
  df <- data.frame(group = grp, Xv)
  fit <- glm(group ~ ., data = df, family = binomial)
  roc_obj <- roc(df$group, predict(fit, type = "response"), quiet = TRUE)
  ci <- ci.auc(roc_obj)
  ggroc(roc_obj, color = "#2C5F8A", linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 1, linetype = 2, color = "grey60", linewidth = 0.3) +
    annotate("text", x = 0.62, y = 0.15, size = 3,
             label = paste0("AUC = ", round(auc(roc_obj), 3),
                            "\n95%CI ", round(ci[1], 2), "-", round(ci[3], 2))) +
    labs(title = title) + theme_pub()
}

hub_genes <- readLines(file.path(OUT_DIR, "04_hub_genes.txt"))

## ============================================================
## 图 1｜Bulk 差异筛选与多队列验证
## ============================================================
p1a <- p1c <- p1d <- NULL
f01 <- file.path(OUT_DIR, "01_GSE134697_DEG_Hippocampus_vs_Neocortex.csv")
if (file.exists(f01)) {
  d01 <- fread(f01)
  p1a <- volcano_plot(d01, "GSE134697: Hippocampus vs Neocortex (paired)")
}
f02a <- file.path(OUT_DIR, "02_GSE190451_DEG_TLE_vs_Control.csv")
f02b <- file.path(OUT_DIR, "02_GSE168375_DEG_Epilepsy_vs_Control.csv")
if (file.exists(f02a)) p1c <- volcano_plot(fread(f02a), "GSE190451: TLE vs Control")
if (file.exists(f02a) && file.exists(f02b)) {
  dA <- fread(f02a); dB <- fread(f02b)
  common <- intersect(dA$gene, dB$gene)
  df_c <- data.frame(x = dA$log2FoldChange[match(common, dA$gene)],
                     y = dB$log2FoldChange[match(common, dB$gene)])
  ct <- cor.test(df_c$x, df_c$y, method = "spearman")
  p1d <- ggplot(df_c, aes(x, y)) +
    geom_point(alpha = 0.25, size = 0.5, color = "grey50") +
    geom_smooth(method = "lm", se = FALSE, color = "#C0392B", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey70", linewidth = 0.3) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 3,
             label = paste0("Spearman rho = ", round(ct$estimate, 3),
                            "\nP = ", format(ct$p.value, digits = 2))) +
    labs(title = "Cross-cohort concordance",
         x = "GSE190451 log2FC (tissue)", y = "GSE168375 log2FC (astrocytes)") +
    theme_pub()
}
panels1 <- list(p1a, p1c, p1d); panels1 <- panels1[!sapply(panels1, is.null)]
if (length(panels1) > 0)
  save_fig(wrap_plots(panels1, nrow = 1) + plot_annotation(tag_levels = "A"),
           "Figure1_bulk_DE_validation", 12, 4.2)

## ============================================================
## 图 2｜机器学习 hub 基因 + ROC
## ============================================================
ml <- fread(file.path(OUT_DIR, "04_ML_selected_genes.csv"))
p2a <- ggplot(ml, aes(x = method, fill = method)) +
  geom_bar(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c(LASSO = "#2C5F8A", SVM_RFE = "#5B8E7D", RF = "#D4A24E")) +
  labs(title = paste0("ML selection (intersection: ", paste(hub_genes, collapse = ", "), ")"),
       x = NULL, y = "Genes selected") + theme_pub()

p2b <- p2c <- p2d <- NULL
try(p2b <- roc_from_csv(file.path(OUT_DIR, "02_GSE168375_vst_normalized.csv"),
                        file.path(OUT_DIR, "02_GSE168375_sample_meta.csv"),
                        hub_genes, "Discovery: GSE168375 (AUC 0.86)"), silent = TRUE)
if (inherits(p2b, "try-error")) p2b <- NULL
try(p2c <- roc_from_csv(file.path(OUT_DIR, "02_GSE190451_normalized.csv"),
                        file.path(OUT_DIR, "02_GSE190451_sample_meta.csv"),
                        hub_genes, "Validation: GSE190451 tissue"), silent = TRUE)
if (inherits(p2c, "try-error")) p2c <- NULL
try(p2d <- roc_from_csv(file.path(OUT_DIR, "08_GSE143272_normalized.csv"),
                        file.path(OUT_DIR, "08_GSE143272_sample_meta.csv"),
                        hub_genes, "Validation: GSE143272 blood"), silent = TRUE)
if (inherits(p2d, "try-error")) p2d <- NULL
panels2 <- list(p2a, p2b, p2c, p2d); panels2 <- panels2[!sapply(panels2, is.null)]
save_fig(wrap_plots(panels2, nrow = 1) + plot_annotation(tag_levels = "A"),
         "Figure2_ML_hub_ROC", 15, 4.2)
## ============================================================
## 图 3｜免疫浸润 + hub 基因与髓系浸润相关性（完整独立版）
## ============================================================

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")
MAIN_DIR <- file.path(BASE_DIR, "figures", "main_figures")
dir.create(MAIN_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
          plot.tag = element_text(face = "bold", size = base_size + 4),
          legend.position = "right",
          strip.background = element_rect(fill = "grey92", color = NA),
          strip.text = element_text(face = "bold"))
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(MAIN_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(MAIN_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat(name, "已保存\n")
}

## 星号标注函数（* P<0.05, ** P<0.01, *** P<0.001）
star_label <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

## ---------- Panel A：免疫浸润 ssGSEA 箱线图 ----------
imm  <- fread(file.path(OUT_DIR, "05_ssGSEA_immune_scores.csv"))
mt1  <- fread(file.path(OUT_DIR, "01_GSE134697_sample_meta.csv"))

imm_l <- melt(imm, id.vars = "cell_type", variable.name = "sample", value.name = "score")
imm_l$sample <- gsub("\\.", "-", imm_l$sample)          # CHKJ.C10 → CHKJ-C10
imm_l$group  <- mt1$tissue[match(imm_l$sample, mt1$sample)]
imm_l <- imm_l[!is.na(group)]
cat("匹配上样本数:", length(unique(imm_l$sample)), "\n")  # 应为 34
stopifnot(nrow(imm_l) > 0)

key_cells <- c("Microglia", "Macrophage", "Monocyte",
               "T_cell_CD4", "T_cell_CD8", "NK_cell")
plot_d <- imm_l[cell_type %in% key_cells]

## 每个细胞类型手动算 Wilcoxon P 值并标注星号
pvals <- plot_d[, {
  w <- wilcox.test(score ~ group)
  .(p = w$p.value,
    y = max(score) + 0.10 * (max(score) - min(score)))
}, by = cell_type]
pvals$label <- star_label(pvals$p)
pvals$x <- 1.5

p3a <- ggplot(plot_d, aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.55) +
  facet_wrap(~cell_type, scales = "free", nrow = 1) +
  scale_fill_manual(values = c(Neocortex = "#7BA7BC", Hippocampus = "#C0392B")) +
  geom_text(data = pvals, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.2) +
  labs(title = "Immune infiltration (ssGSEA, GSE134697)",
       x = NULL, y = "Enrichment score") +
  theme_pub() + theme(legend.position = "none")

## ---------- Panel B：hub 基因与小胶质/巨噬浸润相关性 ----------
fcor <- file.path(OUT_DIR, "05_hub_vs_microglia_macrophage_cor.csv")
p3b <- NULL
if (file.exists(fcor)) {
  cd <- fread(fcor)
  cd_l <- melt(cd, id.vars = "gene",
               measure.vars = c("Microglia_cor", "Macrophage_cor"),
               variable.name = "cell", value.name = "cor")
  p3b <- ggplot(cd_l, aes(gene, cor, fill = cell)) +
    geom_col(position = "dodge", width = 0.65) +
    scale_fill_manual(values = c(Microglia_cor = "#C0392B",
                                 Macrophage_cor = "#7BA7BC"),
                      labels = c("Microglia", "Macrophage")) +
    labs(title = "Hub genes vs myeloid infiltration",
         x = NULL, y = "Spearman r", fill = NULL) +
    theme_pub()
} else cat("未找到 05_hub_vs_microglia_macrophage_cor.csv，跳过 Panel B\n")

## ---------- 拼图并保存 ----------
panels3 <- list(p3a, p3b); panels3 <- panels3[!sapply(panels3, is.null)]
save_fig(wrap_plots(panels3, nrow = 1, widths = c(2, 1.2)) +
           plot_annotation(tag_levels = "A"),
         "Figure3_immune_enrichment", 13, 4.2)

cat("\n输出目录:", MAIN_DIR, "\n")
## ============================================================
## 图 4｜单细胞 LDAM 亚群（需要 06 的 rds）
## ============================================================
fmg <- file.path(OUT_DIR, "06_GSE190452_microglia.rds")
fobj <- file.path(OUT_DIR, "06_GSE190452_annotated_full.rds")
if (file.exists(fmg) && file.exists(fobj)) {
  library(Seurat)
  obj <- readRDS(fobj); mg <- readRDS(fmg)
  mg <- JoinLayers(mg)

  p4a <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE,
                 label.size = 2.6) +
    ggtitle("GSE190452 snRNA-seq") + theme_pub() + theme(legend.position = "none")
  p4b <- ggplot(mg@meta.data, aes(group, LDAM_UCell, fill = group)) +
    geom_violin(scale = "width") + geom_boxplot(width = 0.12, outlier.size = 0.2) +
    scale_fill_manual(values = COL) +
    labs(title = "Microglia LDAM score", x = NULL, y = "UCell score") +
    theme_pub() + theme(legend.position = "none")
  fprop <- file.path(OUT_DIR, "06_LDAm_high_low_proportion_per_sample.csv")
  p4c <- NULL
  if (file.exists(fprop)) {
    pp <- fread(fprop)[LDAM_level == "LDAM-high"]
    p4c <- ggplot(pp, aes(group, proportion, fill = group)) +
      geom_boxplot(width = 0.4) + geom_point(size = 2) +
      scale_fill_manual(values = COL) +
      labs(title = "LDAM-high proportion", x = NULL, y = "Proportion") +
      theme_pub() + theme(legend.position = "none")
  }
  p4d <- FeaturePlot(mg, features = "LDAM_UCell") +
    scale_color_gradient(low = "grey90", high = "#C0392B") +
    ggtitle("LDAM score (microglia)") + theme_pub()
  panels4 <- list(p4a, p4b, p4c, p4d); panels4 <- panels4[!sapply(panels4, is.null)]
  save_fig(wrap_plots(panels4, nrow = 2) + plot_annotation(tag_levels = "A"),
           "Figure4_scRNA_LDAM", 11, 8.5)

  ## ============================================================
  ## 图 5｜LDAM-铁死亡耦联
  ## ============================================================
  p5a <- ggplot(mg@meta.data, aes(LDAM_UCell, Ferroptosis_UCell, color = group)) +
    geom_point(alpha = 0.25, size = 0.4) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
    scale_color_manual(values = COL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 3,
             label = "rho = 0.452, P < 2.2e-16") +
    labs(title = "GSE190452 microglia (discovery)") + theme_pub()
  p5b <- ggplot(mg@meta.data, aes(LDAM_level, Ferroptosis_UCell, fill = LDAM_level)) +
    geom_violin(scale = "width") + geom_boxplot(width = 0.12, outlier.size = 0.2) +
    ggpubr::stat_compare_means(method = "wilcox.test", size = 3) +
    scale_fill_manual(values = COL) +
    labs(title = "Ferroptosis: LDAM-high vs low", x = NULL, y = "Ferroptosis score") +
    theme_pub() + theme(legend.position = "none")

  p5c <- NULL
  fmg7 <- file.path(OUT_DIR, "07_GSE201048_microglia.rds")
  if (file.exists(fmg7)) {
    mg7 <- readRDS(fmg7)
    p5c <- ggplot(mg7@meta.data, aes(LDAM_UCell, Ferroptosis_UCell)) +
      geom_point(alpha = 0.15, size = 0.3, color = "#2C5F8A") +
      geom_smooth(method = "lm", color = "#C0392B", linewidth = 0.7) +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 3,
               label = "rho = 0.317, P < 2.2e-16") +
      labs(title = "GSE201048 microglia (validation, n = 34,479)") + theme_pub()
  }
  panels5 <- list(p5a, p5b, p5c); panels5 <- panels5[!sapply(panels5, is.null)]
  save_fig(wrap_plots(panels5, nrow = 1) + plot_annotation(tag_levels = "A"),
           "Figure5_LDAM_ferroptosis_coupling", 13, 4.2)
} else {
  message("缺少 06 的 rds，跳过图 4、图 5")
}

## ============================================================
## 图 6｜细胞通讯（需要 06 CellChat rds）
## ============================================================
fchat <- file.path(OUT_DIR, "06_CellChat_object.rds")
if (file.exists(fchat) && requireNamespace("CellChat", quietly = TRUE)) {
  library(CellChat)
  chat <- readRDS(fchat)
  pdf(file.path(MAIN_DIR, "Figure6_CellChat.pdf"), width = 13, height = 4.5)
  par(mfrow = c(1, 3))
  netVisual_circle(chat@net$count, vertex.weight = as.numeric(table(chat@idents)),
                   weight.scale = TRUE, label.edge = FALSE, title.name = "Overall network")
  for (pw in c("SPP1", "Cholesterol")) {
    if (pw %in% chat@netP$pathways)
      print(netVisual_aggregate(chat, signaling = pw, layout = "circle"))
  }
  dev.off()
  cat("Figure6_CellChat 已保存\n")
} else message("缺少 CellChat rds，跳过图 6")

## ============================================================
## 图 7｜外周血与动物模型
## ============================================================
p7a <- p7b <- NULL
fblood <- file.path(OUT_DIR, "08_GSE143272_pathway_scores.csv")
fmetaB <- file.path(OUT_DIR, "08_GSE143272_sample_meta.csv")
if (file.exists(fblood) && file.exists(fmetaB)) {
  sb <- fread(fblood); mb <- fread(fmetaB)
  db <- merge(sb, mb, by = "sample")
  db_l <- pivot_longer(db, c(LDAM, Ferroptosis), names_to = "pathway", values_to = "score")
  p7a <- ggplot(db_l, aes(group, score, fill = group)) +
    geom_boxplot(outlier.size = 0.4, width = 0.55) +
    facet_wrap(~pathway, scales = "free") +
    ggpubr::stat_compare_means(method = "wilcox.test", size = 3) +
    scale_fill_manual(values = COL) +
    labs(title = "GSE143272 blood pathway scores", x = NULL, y = "ssGSEA score") +
    theme_pub() + theme(legend.position = "none")
}
fmouse <- file.path(OUT_DIR, "08_GSE88992_pathway_scores.csv")
if (file.exists(fmouse)) {
  sm <- fread(fmouse)
  sm_l <- pivot_longer(sm, c(LDAM, Ferroptosis), names_to = "pathway", values_to = "score")
  p7b <- ggplot(sm_l, aes(group, score, fill = group)) +
    geom_boxplot(outlier.size = 0.4, width = 0.55) + geom_point(size = 1.2) +
    facet_wrap(~pathway, scales = "free") +
    ggpubr::stat_compare_means(method = "wilcox.test", size = 3) +
    scale_fill_manual(values = COL) +
    labs(title = "GSE88992 mouse KA hippocampus", x = NULL, y = "ssGSEA score") +
    theme_pub() + theme(legend.position = "none")
}
panels7 <- list(p7a, p7b); panels7 <- panels7[!sapply(panels7, is.null)]
if (length(panels7) > 0)
  save_fig(wrap_plots(panels7, nrow = 1) + plot_annotation(tag_levels = "A"),
           "Figure7_blood_mouse_validation", 10, 4.2)

cat("\n== 全部主图输出至 ", MAIN_DIR, " ==\n")
print(list.files(MAIN_DIR))
