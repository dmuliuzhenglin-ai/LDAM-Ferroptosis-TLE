# ============================================================================
# Script      : 36_LGALS3_exosome_coupling.R
# Purpose     : Specificity test: LGALS3 co-variation with exosome modules given LDAM score (pseudobulk Spearman partial correlation, with internal reference sets).
# Manuscript  : Supplementary analyses
# Input       : GSE201048 pseudobulk framework (as 34)
# Output      : results/36_LGALS3_exosome_coupling/
# Dependencies: Seurat, dplyr, ggplot2, ppcor, tibble, tidyr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 36_LGALS3_exosome_coupling.R
## 目的：检验 LGALS3 表达是否与外泌体机器模块【特异性】共变
##      LGALS3 ~ 模块评分 | LDAM评分（Spearman 偏相关，ppcor）
##      关键对照：以 LDAM_ref / Ferroptosis_ref 两个参照集为内对照，
##      回答"Gal-3 挂在外泌体臂上"还是"与所有 LDAM 相关通路平行共变"
## 框架：完全复用 34 号脚本（GSE201048 微胶质伪bulk，sample×LDAM分层）
## 输出：results/36_LGALS3_exosome_coupling/
##   36_partial_correlation_modules.csv
##   36_module_scatter.pdf
##   36_summary.txt
## ============================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tibble); library(tidyr)
  library(ppcor); library(ggplot2)
})
set.seed(2024)

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT <- file.path(BASE_DIR, "results/36_LGALS3_exosome_coupling"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## 与 27 号脚本完全一致的基因集（不要改）
gene_sets <- list(
  ESCRT_core        = c("TSG101","VPS4A","VPS4B","VPS25","VPS36","HGS","STAM","STAM2",
                        "CHMP1A","CHMP1B","CHMP2A","CHMP2B","CHMP3","CHMP4A","CHMP4B",
                        "CHMP4C","CHMP5","CHMP6","CHMP7","PDCD6IP"),
  Rab_secretory     = c("RAB27A","RAB27B","RAB11A","RAB11B","RAB35","RAB8A","RAB8B","RAB7A"),
  Ceramide_nSMase   = c("SMPD2","SMPD3","SMPD4"),
  Cargo_loading     = c("YBX1","SYNCRIP","HNRNPA2B1","FMR1","LIN28B","IGF2BP1"),
  SNARE_fusion      = c("SNAP23","STX4","STX6","STX7","VAMP3","VAMP7","VTI1B"),
  Syntenin_ectosome = c("SDCBP","SDCBP2","PDCD6IP","ARRDC1","VPS4A","VPS4B"),
  Exosome_uptake    = c("CD9","CD81","ICAM1","ITGB1","ITGAV","SDC1","SDC4","HSPG2")
)
ref_sets <- list(
  LDAM        = c("PLIN2","DGAT1","DGAT2","ACAT1","ACSL1","CD36","FABP3","HILPDA","SOAT1","LPIN1"),
  Ferroptosis = c("GPX4","ACSL4","LPCAT3","ALOX15","AIFM2","TFRC","FTH1","FTL","SLC7A11","NCOA4")
)

## ---------- 1. 读对象（与 34 一致） ----------
rds <- file.path(BASE_DIR, "results/07_GSE201048_microglia.rds")
stopifnot(file.exists(rds))
mg <- readRDS(rds)
if (inherits(mg[["RNA"]], "Assay5")) mg <- JoinLayers(mg)

## LDAM 评分（与 34 一致：优先复用对象里已有的 LDAM_score）
if (!"LDAM_score" %in% colnames(mg@meta.data)) {
  mx <- as.matrix(GetAssayData(mg, layer = "data"))
  g  <- intersect(ref_sets$LDAM, rownames(mx))
  mg$LDAM_score <- rowMeans(scale(t(mx[g, , drop = FALSE])), na.rm = TRUE)
}
med <- median(mg$LDAM_score, na.rm = TRUE)
mg$LDAM_group <- ifelse(mg$LDAM_score > med, "LDAM_high", "LDAM_low")

## ---------- 2. 伪bulk（sampleName × LDAM_group，与 34 一致） ----------
scol <- "sampleName"
counts <- as.matrix(GetAssayData(mg, layer = "counts"))
pb_key <- paste(mg@meta.data[[scol]], mg$LDAM_group, sep = "__")
pb <- t(rowsum(t(counts), group = pb_key))
keep <- colSums(pb) >= 1e5
pb <- pb[, keep, drop = FALSE]
logcpm <- log2(t(t(pb) / colSums(pb) * 1e6) + 1)
message("伪bulk 数: ", ncol(logcpm))
stopifnot("LGALS3" %in% rownames(logcpm), ncol(logcpm) >= 8)
lgals3  <- logcpm["LGALS3", ]
pb_ldam <- tapply(mg$LDAM_score, pb_key, mean)[colnames(logcpm)]
sapply(mod_list, length)      # 哪一列不是 22
lapply(c(gene_sets, ref_sets), function(g) length(intersect(g, rownames(logcpm))))  # 每模块检出基因数
## ---------- 3. 模块评分（z-score 均值，与 27 号脚本同法） ----------
## 1. 修函数
mod_score <- function(genes, mat){
  g <- intersect(genes, rownames(mat))
  if (length(g) < 3) return(rep(NA_real_, ncol(mat)))
  rowMeans(scale(t(mat[g, , drop = FALSE])), na.rm = TRUE)
}

## 2. 重建（干净版）
mod_list <- lapply(c(gene_sets, ref_sets), mod_score, mat = logcpm)
names(mod_list) <- names(c(gene_sets, ref_sets))
stopifnot(all(sapply(mod_list, length) == ncol(logcpm)))

rownames(mod_mat) <- colnames(logcpm)   # 行 = 22 个伪bulk
colnames(mod_mat) <- names(mod_list)    # 列 = 9 个模块
stopifnot(is.matrix(mod_mat), ncol(mod_mat) == 9, identical(rownames(mod_mat), colnames(logcpm)))
print(mod_mat)
dim(mod_mat)          # 应为 22 9
head(mod_mat[, 1:3])  # 前几行前几列
## ---------- 4. LGALS3 ~ 模块 | LDAM 偏相关 ----------
res <- lapply(colnames(mod_mat), function(md){
  x <- mod_mat[, md]
  if (sd(x, na.rm = TRUE) == 0) return(NULL)
  pc <- tryCatch(pcor.test(lgals3, x, pb_ldam, method = "spearman"), error = function(e) NULL)
  sc <- cor.test(lgals3, x, method = "spearman")
  data.frame(module = md,
             partial_rho  = unname(pc$estimate), partial_P = pc$p.value,
             marginal_rho = unname(sc$estimate), marginal_P = sc$p.value,
             class = ifelse(md %in% names(ref_sets), "reference", "exosome"))
}) %>% bind_rows()
res$FDR <- p.adjust(res$partial_P, "BH")
res <- res %>% arrange(partial_P)
write.csv(res, file.path(OUT, "36_partial_correlation_modules.csv"), row.names = FALSE)
print(res)

## ---------- 5. 散点图 ----------
plotdf <- as.data.frame(mod_mat); plotdf$LGALS3 <- lgals3; plotdf$LDAM <- pb_ldam
plist <- lapply(colnames(mod_mat), function(md){
  r <- res[res$module == md, ]
  ggplot(plotdf, aes(.data[[md]], LGALS3, color = LDAM)) +
    geom_point(size = 2.5, alpha = .85) +
    geom_smooth(method = "lm", se = FALSE, color = "grey40", linetype = 2) +
    scale_color_gradient2(low = "#2E86C8", mid = "grey80", high = "#C0392B", midpoint = med) +
    labs(title = md, subtitle = sprintf("partial rho=%.2f, P=%.3g, FDR=%.3g", r$partial_rho, r$partial_P, r$FDR),
         x = "module score (z)", y = "LGALS3 log2CPM") +
    theme_classic(base_size = 10)
})
pdf(file.path(OUT, "36_module_scatter.pdf"), width = 12, height = 8)
print(patchwork::wrap_plots(plist, ncol = 5)); dev.off()

## ---------- 6. 自动结论 ----------
exo <- res %>% filter(class == "exosome")
ref <- res %>% filter(class == "reference")
exo_sig <- exo %>% filter(FDR < 0.05)
ref_sig <- ref %>% filter(FDR < 0.05)
verdict <- if (nrow(exo_sig) > 0 && nrow(ref_sig) == 0) {
  sprintf("【特异共变】控制 LDAM 后 LGALS3 仍与 %d 个外泌体模块显著共变（%s），而 LDAM/Ferroptosis 参照集自身均不显著 → 支持 Gal-3 特异性挂在外泌体臂上，可写进 Results",
          nrow(exo_sig), paste(exo_sig$module, collapse = ", "))
} else if (nrow(exo_sig) > 0 && nrow(ref_sig) > 0) {
  sprintf("【非特异共变】外泌体模块（%s）与参照集均显著 → Gal-3 与外泌体臂的关联是 LDAM 程序级共变的一部分，写 Discussion 而非 Results", paste(exo_sig$module, collapse = ", "))
} else {
  sprintf("【平行输出】控制 LDAM 后无模块幸存（外泌体最小 P=%.3g）→ Gal-3 与外泌体激活是 LDAM 的平行输出，如实报告，手稿无需改动", min(exo$partial_P, na.rm = TRUE))
}
cat(c("==== 分析 36 结论 ====", verdict, "",
      paste(sprintf("%s (%s): partial rho=%.2f P=%.3g FDR=%.3g",
                    res$module, res$class, res$partial_rho, res$partial_P, res$FDR), collapse = "\n")),
    file = file.path(OUT, "36_summary.txt"), sep = "\n")
cat(readLines(file.path(OUT, "36_summary.txt")), sep = "\n")
