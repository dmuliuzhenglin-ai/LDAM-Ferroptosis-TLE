#!/usr/bin/env Rscript
# ============================================================================
# Script      : 29_miRNA_derepression_panel.R
# Purpose     : Panel test: collective target de-repression for the remaining candidate miRNAs (TargetScan 8.0 conserved targets).
# Manuscript  : Supplementary analyses (miRNA axis)
# Input       : ref/TargetScan_*_targets.csv; GSE168375 logFC
# Output      : results/29_* tables/figures
# Dependencies: tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================
# ===========================================================================
# 29_miRNA_derepression_panel.R
# 收官分析 2：8 个候选 miRNA 的靶基因集体去抑制检验（面板版）
# 背景：脚本 26 已对 miR-34a-5p 完成检验（阴性）。本脚本覆盖其余候选。
# 数据：TargetScan 8.0 保守靶表（已生成，放在 ref/ 目录）：
#   - TargetScan_miR-382-5p_targets.csv   (218 基因)
#   - TargetScan_let-7g-5p_targets.csv    (1207 基因；let-7g/let-7f 同家族同表)
#   - TargetScan_miR-138-5p_targets.csv   (705 基因)
# 注意：miR-191-3p / miR-27b-5p / miR-1912 在 TargetScan 保守预测中不存在
#       （仅有 miR-191-5p / miR-27-3p 臂），如实报告为「无保守预测可检验」
# 检验：GSE168375（31v12）mRNA logFC 分布，靶 vs 非靶单尾 Wilcoxon + 1000 次置换
# 依赖：tidyverse
# ===========================================================================

suppressPackageStartupMessages(library(tidyverse))

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/29_miRNA_derepression", showWarnings = FALSE, recursive = TRUE)

find1 <- function(pattern, exclude = NULL) {
  f <- list.files("data", pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!is.null(exclude)) f <- f[!grepl(exclude, basename(f))]
  if (length(f) == 0) stop("找不到: ", pattern)
  message("使用文件: ", f[1]); f[1]
}

# ---------------------------------------------------------------------------
# 1. GSE168375 表达矩阵 + 逐基因 logFC（与脚本 24/26 完全一致的口径）
# ---------------------------------------------------------------------------
e168 <- read_csv(find1("GSE168375_expr\\.csv$"), show_col_types = FALSE)
smp168 <- setdiff(colnames(e168), "gene")
stopifnot(length(smp168) == 43)
grp168 <- c(rep("TLE", 31), rep("Con", 12))

mat <- as.matrix(e168[, smp168]); rownames(mat) <- e168$gene
# 注意：GSE168375_expr.csv 是 count 矩阵（脚本 24 从 RAW.tar 构建）
# 低表达过滤须在 log2 前
cpm <- t(t(mat) / colSums(mat) * 1e6)
keep <- rowSums(cpm >= 1) >= 10
mat <- log2(mat[keep, ] + 1)
message("检测基因数: ", nrow(mat))

gene_de <- tibble(gene = rownames(mat),
                  logFC = rowMeans(mat[, grp168 == "TLE"]) - rowMeans(mat[, grp168 == "Con"]))

# ---------------------------------------------------------------------------
# 2. 去抑制检验函数（与脚本 26 相同：单尾 Wilcoxon + 1000 置换）
# ---------------------------------------------------------------------------
set.seed(42)
derepress_test <- function(target_file, mir_label) {
  tg <- read_csv(target_file, show_col_types = FALSE)
  tg_genes <- intersect(tg$`Gene Symbol`, gene_de$gene)
  tg_hi <- intersect(tg$`Gene Symbol`[tg$high_confidence], gene_de$gene)
  message(mir_label, ": 靶表 ", nrow(tg), " 基因, 检出 ", length(tg_genes),
          " (高置信 ", length(tg_hi), ")")

  run_one <- function(tset, tag) {
    if (length(tset) < 20) return(tibble(miRNA = mir_label, set = tag, n = length(tset),
                                         median_target = NA, median_nontarget = NA,
                                         P_wilcox = NA, P_perm = NA))
    v_t <- gene_de$logFC[gene_de$gene %in% tset]
    v_n <- gene_de$logFC[!gene_de$gene %in% tset]
    p_w <- wilcox.test(v_t, v_n, alternative = "greater")$p.value
    obs <- median(v_t) - median(v_n)
    perm <- replicate(1000, median(sample(v_n, length(v_t))) - median(v_n))
    p_perm <- mean(perm >= obs)
    tibble(miRNA = mir_label, set = tag, n = length(tset),
           median_target = median(v_t), median_nontarget = median(v_n),
           P_wilcox = p_w, P_perm = p_perm)
  }
  bind_rows(run_one(tg_genes, "all_conserved"), run_one(tg_hi, "high_confidence"))
}

res <- bind_rows(
  derepress_test("ref/TargetScan_miR-382-5p_targets.csv", "miR-382-5p"),
  derepress_test("ref/TargetScan_let-7g-5p_targets.csv", "let-7-5p family (let-7g/let-7f)"),
  derepress_test("ref/TargetScan_miR-138-5p_targets.csv", "miR-138-5p")
)

# 合入 miR-34a-5p 的既有结果（脚本 26）
res26 <- tibble(miRNA = "miR-34a-5p (from script 26)", set = "all_conserved",
                n = 655, median_target = 0.491, median_nontarget = 0.486,
                P_wilcox = 0.31, P_perm = 0.45)
res_all <- bind_rows(res, res26)
res_all$FDR <- p.adjust(res_all$P_wilcox, "BH")
print(res_all, n = Inf)
write_csv(res_all, "results/29_miRNA_derepression/29_derepression_panel.csv")

# 无保守预测的候选（如实记录）
no_ts <- tibble(miRNA = c("miR-191-3p", "miR-27b-5p", "miR-1912"),
                note = "No conserved TargetScan 8.0 prediction exists for this arm; not testable")
write_csv(no_ts, "results/29_miRNA_derepression/29_candidates_without_conserved_predictions.csv")

message("完成。判读：8 候选全部阴性 => miRNA 层面整体为现象级，外泌体路线更坚定；")
message("任一阳性 => 该 miRNA 值得单独跟进荧光素酶实验。")
