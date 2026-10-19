#!/usr/bin/env Rscript
# ============================================================================
# Script      : 28修正.R
# Purpose     : v2 fix of 28_sc_exosome_modules.R (UCell adds _UCell suffix automatically; columns referenced at access time). Supersedes 28.
# Manuscript  : Supplementary analyses (single-cell EV modules)
# Input       : as 28_sc_exosome_modules.R
# Output      : results/28_* tables/figures
# Dependencies: Seurat, UCell, ggplot2, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================
# ===========================================================================
# 28_sc_exosome_modules.R  【v2 修正版】
# 收官分析 1：单细胞层面的外泌体机器模块评分
# 问题：外泌体「摄取/装载臂」激活是否恰好发生在 LDAM-high 小胶质细胞里？
# 队列：GSE201048（3.4 万小胶质细胞）+ GSE190452（如本地有 rds）
# v2 修复：UCell 自动加 _UCell 后缀，不再手动预加；取列时才加
# 依赖：Seurat, UCell, tidyverse, ggplot2
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(tidyverse)
  library(ggplot2)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/28_sc_exosome", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. 基因集（与脚本 27 完全一致；不要加 _UCell 后缀，UCell 自己会加）
# ---------------------------------------------------------------------------
gene_sets <- list(
  ESCRT_core        = c("TSG101","VPS4A","VPS4B","VPS25","VPS36","HGS","STAM","STAM2",
                        "CHMP1A","CHMP1B","CHMP2A","CHMP2B","CHMP3","CHMP4A","CHMP4B",
                        "CHMP4C","CHMP5","CHMP6","CHMP7","PDCD6IP"),
  Rab_secretory     = c("RAB27A","RAB27B","RAB11A","RAB11B","RAB35","RAB8A","RAB8B","RAB7A"),
  Ceramide_nSMase   = c("SMPD2","SMPD3","SMPD4"),
  Cargo_loading     = c("YBX1","SYNCRIP","HNRNPA2B1","FMR1","LIN28B","IGF2BP1"),
  SNARE_fusion      = c("SNAP23","STX4","STX6","STX7","VAMP3","VAMP7","VTI1B"),
  Syntenin_ectosome = c("SDCBP","SDCBP2","PDCD6IP","ARRDC1","VPS4A","VPS4B"),
  Exosome_uptake    = c("CD9","CD81","ICAM1","ITGB1","ITGAV","SDC1","SDC4","HSPG2"),
  LDAM_ref          = c("PLIN2","DGAT1","DGAT2","ACAT1","ACSL1","CD36","FABP3","HILPDA","SOAT1","LPIN1"),
  Ferroptosis_ref   = c("GPX4","ACSL4","LPCAT3","ALOX15","AIFM2","TFRC","FTH1","FTL","SLC7A11","NCOA4")
)

# ---------------------------------------------------------------------------
# 2. 通用分析函数
# ---------------------------------------------------------------------------
analyze_obj <- function(obj, dataset_label) {
  message("===== ", dataset_label, " =====")
  message("细胞数: ", ncol(obj))

  # 2a. UCell 打分（MulticoreParam 警告在 Windows 上可忽略，自动单线程）
  obj <- AddModuleScore_UCell(obj, features = gene_sets, ncores = 1)

  meta <- obj@meta.data
  score_cols <- paste0(names(gene_sets), "_UCell")   # UCell 自动加的后缀
  missing <- setdiff(score_cols, colnames(meta))
  if (length(missing) > 0) stop("UCell 列缺失: ", paste(missing, collapse = ", "))

  # 2b. LDAM-high/low 分层：优先用对象里已有的 LDAM 评分，否则用 LDAM_ref
  ldam_col <- grep("LDAM|LDam", colnames(meta), value = TRUE, ignore.case = TRUE)
  ldam_col <- ldam_col[!ldam_col %in% c("LDAM_ref_UCell")][1]
  if (is.na(ldam_col)) {
    message("未找到已有 LDAM 评分列，用 LDAM_ref_UCell 分层")
    ldam_col <- "LDAM_ref_UCell"
  } else message("沿用已有 LDAM 评分列: ", ldam_col)

  meta$LDAM_group <- ifelse(meta[[ldam_col]] > median(meta[[ldam_col]], na.rm = TRUE),
                            "LDAM_high", "LDAM_low")
  print(table(meta$LDAM_group))

  # 2c. 组间 Wilcoxon 比较（除分层列自身外）
  test_cols <- setdiff(score_cols, ldam_col)
  res <- bind_rows(lapply(test_cols, function(cc) {
    x <- meta[[cc]][meta$LDAM_group == "LDAM_high"]
    y <- meta[[cc]][meta$LDAM_group == "LDAM_low"]
    tibble(dataset = dataset_label, module = str_remove(cc, "_UCell"),
           median_high = median(x, na.rm = TRUE), median_low = median(y, na.rm = TRUE),
           delta = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
           P = tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_))
  }))
  res$FDR <- p.adjust(res$P, "BH")
  print(res, n = Inf)
  list(res = res, meta = meta, ldam_col = ldam_col)
}

# ---------------------------------------------------------------------------
# 3. GSE201048 小胶质细胞
# ---------------------------------------------------------------------------
f048 <- list.files(".", pattern = "07_GSE201048_microglia\\.rds$",
                   recursive = TRUE, full.names = TRUE)[1]
stopifnot(!is.na(f048))
message("使用对象: ", f048)
mg048 <- readRDS(f048)
r048 <- analyze_obj(mg048, "GSE201048_microglia")

# ---------------------------------------------------------------------------
# 4. GSE190452（TLE-HS 单核，如本地已有处理好的 Seurat rds 则一并跑）
# ---------------------------------------------------------------------------
f452 <- list.files(".", pattern = "GSE190452.*\\.rds$", recursive = TRUE, full.names = TRUE)[1]
r452 <- NULL
if (!is.na(f452)) {
  message("找到 GSE190452 对象: ", f452)
  obj452 <- readRDS(f452)
  ann_cols <- c("celltype", "cell_type", "annotation", "cluster")
  ann_col <- intersect(ann_cols, colnames(obj452@meta.data))[1]
  if (!is.na(ann_col) && any(grepl("micro|MG", obj452@meta.data[[ann_col]], ignore.case = TRUE))) {
    mg452 <- subset(obj452, cells = rownames(obj452@meta.data)[grepl("micro|MG", obj452@meta.data[[ann_col]], ignore.case = TRUE)])
    r452 <- analyze_obj(mg452, "GSE190452_microglia")
  } else {
    message("未识别小胶质细胞注释列，请把 obj452 的小胶质细胞子集存成 rds 后告知路径")
  }
} else {
  message("未找到 GSE190452 rds，跳过（可后补）")
}

# ---------------------------------------------------------------------------
# 5. 汇总输出
# ---------------------------------------------------------------------------
all_res <- bind_rows(r048$res, if (!is.null(r452)) r452$res)
write_csv(all_res, "results/28_sc_exosome/28_sc_exosome_module_tests.csv")

plot_meta <- r048$meta %>%
  select(LDAM_group, all_of(paste0(names(gene_sets), "_UCell"))) %>%
  pivot_longer(-LDAM_group, names_to = "module", values_to = "score") %>%
  mutate(module = str_remove(module, "_UCell"))

p <- ggplot(plot_meta, aes(LDAM_group, score, fill = LDAM_group)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.size = 0.2) +
  facet_wrap(~ module, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("LDAM_high" = "#C9828B", "LDAM_low" = "#7A9CC4")) +
  labs(title = "Exosome machinery module scores in LDAM-high vs LDAM-low microglia (GSE201048)",
       x = NULL, y = "UCell score") +
  theme(legend.position = "none", strip.text = element_text(size = 8))
ggsave("results/28_sc_exosome/28_sc_exosome_violin.pdf", p, width = 11, height = 5.5)

message("完成。判读要点：")
message("- 若 Exosome_uptake / Cargo_loading 在 LDAM-high 显著更高(FDR<0.05)，")
message("  第四部分与主线焊死：机器激活恰好发生在 LDAM-high 小胶质细胞里")
message("- Ferroptosis_ref 应显著更高（主文已知结论，作阳性对照）")
