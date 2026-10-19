# ============================================================================
# Script      : 02b_GSE190451_fix.R
# Purpose     : Re-run of the GSE190451 part (3 TLE vs 3 control) plus log2FC direction concordance with GSE168375.
# Manuscript  : Methods 2.2 / Results 3.1
# Input       : GSE190451_TLE_samples.txt.gz; results of 02
# Output      : results/02b_* tables and figures
# Dependencies: DESeq2, data.table, ggplot2, limma
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 02b_GSE190451_fix.R
## 补跑 Part A：GSE190451（3 TLE vs 3 对照，颞叶皮层 bulk）
## + 与 GSE168375 的 log2FC 方向一致性检查
## ============================================================
library(data.table)
library(limma)
library(ggplot2)

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
res_dir  <- file.path(BASE_DIR, "results")
fig_dir  <- file.path(BASE_DIR, "figures")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- 1. 定位并读取 GSE190451 文件 ----------
f451 <- list.files(BASE_DIR, pattern = "GSE190451.*\\.(txt|tsv|csv)(\\.gz)?$",
                   recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
cat("找到文件：\n"); print(f451)
stopifnot(length(f451) >= 1)

d451 <- fread(f451[1])
cat("列名：\n"); print(names(d451))

sample_cols <- grep("^(Normal|TLE)_", names(d451), value = TRUE)
stopifnot(length(sample_cols) == 6)
cat("样本列：", paste(sample_cols, collapse = ", "), "\n")

## 基因名用 Gene_name 列；去掉空名和重复（取各样本之和合并）
d451 <- d451[!is.na(Gene_name) & Gene_name != ""]
mat <- as.matrix(d451[, ..sample_cols])
storage.mode(mat) <- "numeric"
mat <- rowsum(mat, group = d451$Gene_name, reorder = FALSE)
cat("基因数：", nrow(mat), "；样本数：", ncol(mat), "\n")

group <- factor(ifelse(grepl("^Normal", sample_cols), "Control", "TLE"),
                levels = c("Control", "TLE"))
metaA <- data.frame(sample = sample_cols, group = group)
print(table(metaA$group))

## ---------- 2. 差异分析（FPKM 型小数值 → log2 + limma） ----------
cat("数值范围：min =", min(mat), "，max =", max(mat), "\n")

if (all(mat == round(mat))) {
  ## 万一其实是整数计数，则用 DESeq2
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(round(mat), metaA, ~ group)
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds, contrast = c("group", "TLE", "Control")))
  resA <- data.table(gene = rownames(res), log2FoldChange = res$log2FoldChange,
                     pvalue = res$pvalue, padj = res$padj)
} else {
  expr <- log2(mat + 0.5)
  fit <- lmFit(expr, model.matrix(~ group))
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "groupTLE", number = Inf, sort.by = "none")
  resA <- data.table(gene = rownames(tt), log2FoldChange = tt$logFC,
                     pvalue = tt$P.Value, padj = tt$adj.P.Val)
}

resA[, direction := fifelse(padj < 0.05 & log2FoldChange >  0.585, "Up",
                            fifelse(padj < 0.05 & log2FoldChange < -0.585, "Down", "NS"))]
setorder(resA, padj)
fwrite(resA, file.path(res_dir, "02_GSE190451_DEG_TLE_vs_Control.csv"))
fwrite(metaA, file.path(res_dir, "02_GSE190451_sample_meta.csv"))
cat("== GSE190451 DEG 数 ==\n"); print(table(resA$direction))
cat("== 前 15 个最显著基因 ==\n"); print(head(resA, 15))

## ---------- 3. 火山图 ----------
pdf(file.path(fig_dir, "02_volcano_GSE190451.pdf"), width = 7, height = 6)
print(ggplot(resA, aes(log2FoldChange, -log10(pvalue), color = direction)) +
        geom_point(alpha = 0.6, size = 1.2) +
        scale_color_manual(values = c(Up = "#C0392B", Down = "#2C5F8A", NS = "grey80")) +
        geom_vline(xintercept = c(-0.585, 0.585), linetype = 2, color = "grey50") +
        geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
        theme_classic() + labs(title = "GSE190451: TLE vs Control (temporal cortex)"))
dev.off()

## ---------- 4. 与 GSE168375 的 log2FC 一致性检查 ----------
resB <- fread(file.path(res_dir, "02_GSE168375_DEG_Epilepsy_vs_Control.csv"))
common <- intersect(resA$gene, resB$gene)
cat("\n两队列重叠基因数：", length(common), "\n")
mA <- resA[match(common, gene)]$log2FoldChange
mB <- resB[match(common, gene)]$log2FoldChange
ct <- cor.test(mA, mB, method = "spearman")
cat("== 方向一致性（Spearman）==\n")
cat("rho =", round(ct$estimate, 3), "，P =", format(ct$p.value, digits = 3), "\n")

pdf(file.path(fig_dir, "02_log2FC_consistency_GSE190451_vs_GSE168375.pdf"), width = 6, height = 6)
plot(mA, mB, pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.3),
     xlab = "GSE190451 log2FC (tissue)", ylab = "GSE168375 log2FC (astrocytes)",
     main = paste0("Spearman rho = ", round(ct$estimate, 3)))
abline(h = 0, v = 0, col = "grey60", lty = 2)
abline(lm(mB ~ mA), col = "#C0392B")
dev.off()

cat("\n全部完成！\n")
