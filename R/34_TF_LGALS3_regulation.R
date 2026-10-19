# ============================================================================
# Script      : 34_TF_LGALS3_regulation.R
# Purpose     : Partial correlation of TF activity with LGALS3 expression controlling LDAM score (focus: PPARG / NR1H2 / NR1H3).
# Manuscript  : Supplementary analyses (upstream regulators)
# Input       : GSE201048 Seurat object (as 30); TF activities from 30 environment
# Output      : results/34_TF_LGALS3_regulation/
# Dependencies: OmnipathR, Seurat, decoupleR, dplyr, ggplot2, ppcor, tibble, tidyr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================

## ============================================================
## 34_TF_LGALS3_regulation.R
## 目的：检验「TF 活性 → LGALS3 表达」在控制 LDAM 评分后是否仍成立
##      （partial correlation，ppcor），重点：PPARG / NR1H2 / NR1H3
## 依赖：30_TF_activity_LDAM.R 已跑通的环境（Seurat v5 + decoupleR + dorothea_hs）
## 输入：GSE201048 Seurat 对象（与 30 号脚本同一路径）
## 输出：results/34_TF_LGALS3_regulation/
##   34_partial_correlation_all_TF.csv   全 TF 偏相关排序表
##   34_PPARG_focus_scatter.pdf          重点 TF 散点图
##   34_summary.txt                      自动结论
## ============================================================
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# install.packages(c("ppcor","patchwork"))

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tibble); library(tidyr)
  library(decoupleR); library(OmnipathR); library(ppcor); library(ggplot2)
})
set.seed(2024)

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT <- file.path(BASE_DIR, "results/34_TF_LGALS3_regulation"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

LDAM_ref <- c("PLIN2","DGAT2","HILPDA","SOAT1","FITM2","LPL","CD36","APOE","SPP1","GPNMB")

## ---------- 1. 读对象（与 28/30 号脚本一致） ----------
rds <- file.path(BASE_DIR, "results/07_GSE201048_microglia.rds")
stopifnot(file.exists(rds))
message("加载: ", rds)

if (inherits(mg[["RNA"]], "Assay5")) mg <- JoinLayers(mg)   # v5 多层才需要合并；v4 跳过

## LDAM 分层（与 28/30 一致：UCell 若已存在直接用，否则 z-score 均值）
if (!"LDAM_ref_UCell" %in% colnames(mg@meta.data)) {
  mx <- as.matrix(GetAssayData(mg, layer = "data"))
  g  <- intersect(LDAM_ref, rownames(mx))
  z  <- scale(t(mx[g, , drop = FALSE]))
  mg$LDAM_score <- rowMeans(z, na.rm = TRUE)
} else {
  mg$LDAM_score <- mg$LDAM_ref_UCell
}
med <- median(mg$LDAM_score, na.rm = TRUE)
mg$LDAM_group <- ifelse(mg$LDAM_score > med, "LDAM_high", "LDAM_low")
stopifnot(all(c("LDAM_high","LDAM_low") %in% mg$LDAM_group))

## ---------- 2. 伪bulk（sample × LDAM_group），与 30 号脚本一致 ----------
scol <- "sampleName"   # GSE201048 的患者列（11 例）
counts <- as.matrix(GetAssayData(mg, layer = "counts"))
pb_key <- paste(mg@meta.data[[scol]], mg$LDAM_group, sep = "__")
pb <- rowsum(t(counts), group = pb_key)          # pseudobulk × gene
pb <- t(pb)
keep <- colSums(pb) >= 1e5                        # 每个伪bulk的总UMI（按列求和）
pb <- pb[, keep, drop = FALSE]
logcpm <- log2(t(t(pb) / colSums(pb) * 1e6) + 1)
message("伪bulk 数: ", ncol(logcpm))
stopifnot("LGALS3" %in% rownames(logcpm), ncol(logcpm) >= 10)

lgals3 <- logcpm["LGALS3", ]
# 每个伪bulk的 LDAM 连续评分 = 其成员细胞 LDAM_score 均值
pb_ldam <- tapply(mg$LDAM_score, pb_key, mean)[colnames(logcpm)]
stopifnot(!any(is.na(pb_ldam)))

## ---------- 3b. 整理 TF 活性矩阵（防御性写法） ----------
message("acts 列名: ", paste(colnames(acts), collapse = ", "))
print(table(acts$statistic))

## 统一列名：有的版本 statistic 值是 "ulm"，有的直接只有一行 score
acts2 <- as.data.frame(acts)
if ("statistic" %in% colnames(acts2) && "ulm" %in% acts2$statistic) {
  acts2 <- acts2[acts2$statistic == "ulm", ]
}
stopifnot(all(c("source","condition","score") %in% colnames(acts2)))

tf_wide <- acts2 %>%
  dplyr::select(condition, source, score) %>%
  dplyr::distinct(condition, source, .keep_all = TRUE) %>%
  pivot_wider(names_from = source, values_from = score)

tf_mat <- as.matrix(column_to_rownames(as.data.frame(tf_wide), "condition"))
tf_mat <- tf_mat[colnames(logcpm), , drop = FALSE]
stopifnot(identical(rownames(tf_mat), colnames(logcpm)))
message("TF 矩阵: ", nrow(tf_mat), " 伪bulk × ", ncol(tf_mat), " TF")
## ---------- 4. Partial correlation：LGALS3 ~ TF | LDAM ----------
res <- lapply(colnames(tf_mat), function(tf){
  x <- tf_mat[, tf]
  if (sd(x, na.rm = TRUE) == 0) return(NULL)
  pc <- tryCatch(pcor.test(lgals3, x, pb_ldam, method = "spearman"), error = function(e) NULL)
  sc <- cor.test(lgals3, x, method = "spearman")
  data.frame(TF = tf,
             partial_rho   = unname(pc$estimate),
             partial_P     = pc$p.value,
             marginal_rho  = unname(sc$estimate),
             marginal_P    = sc$p.value)
}) %>% bind_rows()
res$FDR <- p.adjust(res$partial_P, "BH")
res <- res %>% arrange(partial_P)
res$focus <- res$TF %in% c("PPARG","NR1H2","NR1H3","CEBPA","SREBF1","RELA","NFKB1")
write.csv(res, file.path(OUT, "34_partial_correlation_all_TF.csv"), row.names = FALSE)

## ---------- 5. 重点 TF 散点图 ----------
focus_tf <- intersect(c("PPARG","NR1H2","NR1H3","CEBPA","SREBF1"), colnames(tf_mat))
plotdf <- as.data.frame(tf_mat[, focus_tf, drop = FALSE])
plotdf$LGALS3 <- lgals3; plotdf$LDAM <- pb_ldam; plotdf$pb <- rownames(plotdf)
plist <- lapply(focus_tf, function(tf){
  r <- res[res$TF == tf, ]
  ggplot(plotdf, aes(.data[[tf]], LGALS3, color = LDAM)) +
    geom_point(size = 2.5, alpha = .85) +
    geom_smooth(method = "lm", se = FALSE, color = "grey40", linetype = 2) +
    scale_color_gradient2(low = "#2E86C8", mid = "grey80", high = "#C0392B", midpoint = med) +
    labs(title = tf,
         subtitle = sprintf("partial rho=%.2f, P=%.3g, FDR=%.3g", r$partial_rho, r$partial_P, r$FDR),
         x = "TF activity (ULM)", y = "LGALS3 log2CPM") +
    theme_classic(base_size = 11)
})
pdf(file.path(OUT, "34_PPARG_focus_scatter.pdf"), width = 10, height = 7)
print(patchwork::wrap_plots(plist, ncol = 3)); dev.off()

{
  pparg <- res[res$TF == "PPARG", ]
  verdict <- if (nrow(pparg) == 0) {
    "PPARG 不在 DoRothEA A-C 网络中"
  } else if (pparg$FDR < 0.05) {
    sprintf("【强支持】PPARG 活性在控制 LDAM 后仍显著预测 LGALS3（partial rho=%.2f, FDR=%.3g）→ 支持 PPARG→Gal-3 轴，实验5设计有分子级依据", pparg$partial_rho, pparg$FDR)
  } else if (pparg$partial_P < 0.05) {
    sprintf("【弱支持】PPARG 名义显著但未过 FDR（partial rho=%.2f, P=%.3g, FDR=%.3g）→ 只能写 nominal，实验5依据降级为探索性", pparg$partial_rho, pparg$partial_P, pparg$FDR)
  } else {
    sprintf("【阴性】控制 LDAM 后 PPARG 不再预测 LGALS3（partial rho=%.2f, P=%.3g）→ PPARG 与 LGALS3 的相关由 LDAM 状态解释，需回到 NR1H2/3 或 CEBPA 候选", pparg$partial_rho, pparg$partial_P)
  }

  top5 <- res %>% slice_head(n = 5) %>%
    mutate(line = sprintf("%s: partial rho=%.2f P=%.3g FDR=%.3g", TF, partial_rho, partial_P, FDR)) %>% pull(line)
  cat(c("==== 分析 34 结论 ====", verdict, "", "偏相关 Top5:", top5,
        "", "注：n = 伪bulk数（22），统计效能有限；FDR<0.05 方可写进 Results"),
      file = file.path(OUT, "34_summary.txt"), sep = "\n")
  cat(readLines(file.path(OUT, "34_summary.txt")), sep = "\n")
}
