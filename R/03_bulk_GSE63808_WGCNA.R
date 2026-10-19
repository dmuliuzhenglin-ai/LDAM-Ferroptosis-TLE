# ============================================================================
# Script      : 03_bulk_GSE63808_WGCNA.R
# Purpose     : WGCNA co-expression modules in GSE63808 (129 TLE hippocampi) and correlation with LDAM/ferroptosis ssGSEA scores.
# Manuscript  : Methods 2.2 / Results 3.3
# Input       : GSE63808 series matrix; genesets/ from 00
# Output      : results/03_* module tables; module-trait figures
# Dependencies: GEOquery, GSVA, WGCNA, data.table, ggplot2, limma, stringr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 03_bulk_GSE63808_WGCNA.R
## GSE63808（129 例 TLE 海马，无对照，HumanHT-12 V3 芯片）
## WGCNA 共表达模块分析 + LDAM/铁死亡 ssGSEA 评分关联
## 修正版：本地读文件 + 强制单线程（修复 Windows 并行报错）
## ============================================================

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results"); FIG_DIR <- file.path(BASE_DIR, "figures")
GS_DIR   <- file.path(BASE_DIR, "genesets")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(GEOquery); library(data.table); library(stringr); library(limma)
  library(WGCNA); library(GSVA); library(ggplot2)
})
disableWGCNAThreads()          # 【关键】关闭多线程，防止 Windows 并行报错
options(stringsAsFactors = FALSE)

## ---------- 0. 读基因集 ----------
LDAM_genes        <- readLines(file.path(GS_DIR, "LDAM_geneset.txt"))
Ferroptosis_genes <- readLines(file.path(GS_DIR, "Ferroptosis_geneset.txt"))
cat("LDAM 基因数：", length(LDAM_genes), "；铁死亡基因数：", length(Ferroptosis_genes), "\n")

## ---------- 1. 读取本地 series matrix ----------
f63808 <- file.path(BASE_DIR, "GSE63808", "GSE63808_series_matrix.txt.gz")
stopifnot(file.exists(f63808))
gse  <- getGEO(filename = f63808)      # 读本地文件时直接返回 ExpressionSet
eset <- gse
expr <- exprs(eset)
cat("表达矩阵：", nrow(expr), "探针 x", ncol(expr), "样本\n")

## 如最大值异常大则说明未取 log2，做一次 log2 变换
if (max(expr, na.rm = TRUE) > 100) {
  expr <- log2(expr + 1)
  cat("已执行 log2 变换\n")
}

## ---------- 2. 探针 → 基因 symbol ----------
fd <- fData(eset)
sym_col <- grep("symbol|gene.?symbol|gene_assignment", names(fd),
                ignore.case = TRUE, value = TRUE)[1]

if (!is.na(sym_col)) {
  cat("使用 series matrix 自带的注释列：", sym_col, "\n")
  symbols <- as.character(fd[[sym_col]])
} else {
  ## 备用：联网/本地读 GPL6947 注释
  cat("series matrix 无 symbol 列，尝试读取 GPL6947 注释……\n")
  gpl_local <- file.path(BASE_DIR, "GSE63808", "GPL6947.annot.gz")
  gpl <- tryCatch(
    if (file.exists(gpl_local)) getGEO(filename = gpl_local) else getGEO("GPL6947", AnnotGPL = TRUE),
    error = function(e) NULL
  )
  if (is.null(gpl)) {
    stop("GPL6947 注释下载失败（网络问题）。请用浏览器打开：\n",
         "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6947/annot/GPL6947.annot.gz\n",
         "下载后放到：", gpl_local, " ，然后重新 Source 本脚本。")
  }
  gpl_tb <- Table(gpl)
  id_col  <- grep("^ID$", names(gpl_tb), value = TRUE)[1]
  gs_col  <- grep("symbol", names(gpl_tb), ignore.case = TRUE, value = TRUE)[1]
  symbols <- gpl_tb[[gs_col]][match(rownames(expr), gpl_tb[[id_col]])]
}

symbols <- str_split_fixed(symbols, " /// ", 2)[, 1]   # 多映射探针取第一个
symbols <- str_trim(symbols)
keep <- !is.na(symbols) & symbols != ""
expr <- expr[keep, ]; symbols <- symbols[keep]
cat("成功注释的探针：", length(symbols), "/", nrow(fd), "\n")

## 同基因多探针合并（取均值），得到基因 x 样本矩阵
expr_gene <- avereps(expr, ID = symbols)
cat("基因水平矩阵：", nrow(expr_gene), "基因 x", ncol(expr_gene), "样本\n")

## ---------- 3. ssGSEA：LDAM / 铁死亡评分（作为 WGCNA 性状） ----------
gsets <- list(LDAM = intersect(LDAM_genes, rownames(expr_gene)),
              Ferroptosis = intersect(Ferroptosis_genes, rownames(expr_gene)))
cat("进入打分的基因：LDAM", length(gsets$LDAM), "；Ferroptosis", length(gsets$Ferroptosis), "\n")

gsva_par <- gsvaParam(exprData = as.matrix(expr_gene), geneSets = gsets)
ss <- gsva(gsva_par, verbose = FALSE)          # 2 x 129
fwrite(as.data.table(t(ss), keep.rownames = "sample"),
       file.path(OUT_DIR, "03_GSE63808_ssGSEA_scores.csv"))

## ---------- 4. 筛选高变异基因构建共表达网络 ----------
mad_v <- apply(expr_gene, 1, mad)
top_n <- min(8000, nrow(expr_gene))
datExpr <- t(expr_gene[order(mad_v, decreasing = TRUE)[1:top_n], ])
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
datTraits <- data.frame(LDAM = ss["LDAM", rownames(datExpr)],
                        Ferroptosis = ss["Ferroptosis", rownames(datExpr)])
stopifnot(!any(is.na(datTraits)))

## ---------- 5. 软阈值 ----------
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,
                         networkType = "signed")
pdf(file.path(FIG_DIR, "03_soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Power", ylab = "Scale-free fit (signed R^2)", type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "#C0392B"); abline(h = 0.8, col = "grey50", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Power",
     ylab = "Mean connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "#2C5F8A")
dev.off()

power <- sft$powerEstimate
if (is.na(power) || length(power) == 0) {
  power <- 6
  cat("！！ 自动软阈值未达到 R2>0.8，采用默认 power = 6（signed 网络常用值）\n")
}
cat("== 选定软阈值 power =", power, "==\n")

## ---------- 6. 构建网络 + 模块划分 ----------
net <- blockwiseModules(datExpr, power = power, networkType = "signed",
                        TOMType = "signed", minModuleSize = 30,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        saveTOMs = FALSE, verbose = 3)
mod_colors <- labels2colors(net$colors)
mod_tab <- table(mod_colors)
cat("== 模块数量（不含 grey）==\n"); print(mod_tab)

## 模块-性状相关
MEs <- moduleEigengenes(datExpr, colors = mod_colors)$eigengenes
MEs <- orderMEs(MEs)
mt_cor  <- cor(MEs, datTraits, use = "p")
mt_p    <- corPvalueStudent(mt_cor, nrow(datExpr))
mt_out <- data.frame(module = rownames(mt_cor),
                     LDAM_r = mt_cor[, "LDAM"], LDAM_P = mt_p[, "LDAM"],
                     Ferro_r = mt_cor[, "Ferroptosis"], Ferro_P = mt_p[, "Ferroptosis"])
mt_out <- mt_out[order(-abs(mt_out$LDAM_r)), ]
print(mt_out)
fwrite(mt_out, file.path(OUT_DIR, "03_WGCNA_module_trait_correlation.csv"))

pdf(file.path(FIG_DIR, "03_module_trait_heatmap.pdf"), width = 7, height = 9)
labeledHeatmap(Matrix = mt_cor, xLabels = names(datTraits),
               yLabels = names(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50), setStdMargins = FALSE,
               textMatrix = paste0(signif(mt_cor, 2), "\n(",
                                   signif(mt_p, 1), ")"),
               main = "Module-Trait relationships")
dev.off()

## ---------- 7. 锁定 LDAM 最高相关模块，导出 hub 候选 ----------
hub_mod <- gsub("^ME", "", mt_out$module[1])
cat("== LDAM 最高相关模块：", hub_mod,
    "（r =", round(mt_out$LDAM_r[1], 3),
    "，P =", format(mt_out$LDAM_P[1], digits = 3), "）==\n")

in_mod <- mod_colors == hub_mod
mm <- cor(datExpr[, in_mod], MEs[[paste0("ME", hub_mod)]], use = "p")
gs_ldam <- cor(datExpr[, in_mod], datTraits$LDAM, use = "p")
hub_df <- data.table(gene = colnames(datExpr)[in_mod],
                     module = hub_mod,
                     module_membership = as.numeric(mm),
                     LDAM_gene_significance = as.numeric(gs_ldam))
hub_df <- hub_df[order(-abs(hub_df$module_membership)), ]
fwrite(hub_df, file.path(OUT_DIR, "03_WGCNA_LDAM_hub_module_genes.csv"))
writeLines(hub_df$gene, file.path(OUT_DIR, "03_WGCNA_LDAM_hub_module_genes.txt"))
setorder(hub_df, -abs(module_membership))
fwrite(hub_df, file.path(OUT_DIR, "03_WGCNA_LDAM_hub_module_genes.csv"))
writeLines(hub_df$gene, file.path(OUT_DIR, "03_WGCNA_LDAM_hub_module_genes.txt"))

pdf(file.path(FIG_DIR, "03_MM_vs_GS_scatter.pdf"), width = 6, height = 6)
plot(as.numeric(mm), as.numeric(gs_ldam), pch = 16, cex = 0.6,
     col = rgb(0.2, 0.4, 0.6, 0.5),
     xlab = paste0("Module Membership (", hub_mod, ")"),
     ylab = "Gene Significance for LDAM score",
     main = paste0("Hub module: ", hub_mod))
abline(lm(as.numeric(gs_ldam) ~ as.numeric(mm)), col = "#C0392B")
dev.off()

cat("\n== 03 完成！模块-性状相关见 03_WGCNA_module_trait_correlation.csv ==\n")
