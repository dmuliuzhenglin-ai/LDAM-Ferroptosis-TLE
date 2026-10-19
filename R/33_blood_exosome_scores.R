#!/usr/bin/env Rscript
# ============================================================================
# Script      : 33_blood_exosome_scores.R
# Purpose     : Detectability of exosome-machinery modules in peripheral blood (GSE143272); per-gene blood-accessible marker candidates.
# Manuscript  : Supplementary analyses (peripheral blood)
# Input       : GSE143272 normalized matrix (from 08)
# Output      : results/33_* tables/figures
# Dependencies: ggplot2, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================
# ===========================================================================
# 33_blood_exosome_scores.R
# 收官分析 3：外泌体机器模块在外周血的可检测性
# 数据：GSE143272（34 未用药 / 57 已用药 / 51 健康对照，Illumina BeadArray）
# 预期：主文已报 bulk 血液 LDAM/铁死亡评分阴性，本分析大概率同样阴性
#       ——如实报告；同时单列外泌体 hub 基因（CD9 等）的逐基因结果，
#       那才是血液可及标志物的候选
# 依赖：tidyverse, ggplot2
# ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/33_blood_exosome", showWarnings = FALSE, recursive = TRUE)

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

module_scores <- function(mat, sets) {
  bind_rows(lapply(names(sets), function(nm) {
    g <- intersect(sets[[nm]], rownames(mat))
    if (length(g) < 3) return(tibble(sample = colnames(mat), module = nm, score = NA_real_, n_genes = length(g)))
    z <- t(scale(t(mat[g, , drop = FALSE]))); z[is.na(z)] <- 0
    tibble(sample = colnames(mat), module = nm, score = colMeans(z), n_genes = length(g))
  }))
}

# ---------------------------------------------------------------------------
# 1. 读血液归一化矩阵
# ---------------------------------------------------------------------------
f_e <- list.files("results", pattern = "08_GSE143272_normalized\\.csv$",
                  recursive = TRUE, full.names = TRUE)[1]
stopifnot(!is.na(f_e))
message("使用文件: ", f_e)
e <- read_csv(f_e, show_col_types = FALSE)
sym_col <- names(e)[1]
message("首列: ", sym_col, " | 样本数: ", ncol(e) - 1)

mat <- as.matrix(e[, -1]); rownames(mat) <- e[[sym_col]]
# 若存在重复基因名，取均值
if (anyDuplicated(rownames(mat))) {
  message("合并重复基因名（取均值）")
  mat <- as.data.frame(mat) %>% rownames_to_column("g") %>%
    group_by(g) %>% summarise(across(everything(), mean), .groups = "drop") %>%
    column_to_rownames("g") %>% as.matrix()
}

# ---------------------------------------------------------------------------
# 2. 分组（34 drug-naive / 57 treated / 51 healthy）
# ---------------------------------------------------------------------------
f_m <- list.files("results", pattern = "08_GSE143272.*(meta|sample|group|pheno)",
                  recursive = TRUE, full.names = TRUE)[1]
stopifnot(!is.na(f_m))
message("使用分组文件: ", f_m)
meta <- read_csv(f_m, show_col_types = FALSE)
glimpse(meta)

# 自动识别样本列与分组列
scol <- names(meta)[sapply(meta, function(x) sum(as.character(x) %in% colnames(mat)) > 30)][1]
gcol <- names(meta)[sapply(meta, function(x)
  any(grepl("naive|treated|healthy|control", as.character(x), ignore.case = TRUE)))][1]
message("样本列: ", scol, " | 分组列: ", gcol)
stopifnot(!is.na(scol), !is.na(gcol))

grp <- set_names(as.character(meta[[gcol]]), as.character(meta[[scol]]))
print(table(grp))
# 统一标签（本 meta 仅含 34 未用药 + 51 健康对照）
grp <- set_names(ifelse(grepl("Healthy", as.character(meta[[gcol]]), ignore.case = TRUE),
                        "Healthy", "Naive"),
                 as.character(meta[[scol]]))
print(table(grp))
stopifnot(sum(grp == "Naive") == 34, sum(grp == "Healthy") == 51)

# ---------------------------------------------------------------------------
# 3. 模块评分 + Naive vs Healthy
# ---------------------------------------------------------------------------
scores <- module_scores(mat, gene_sets)
scores$grp <- grp[scores$sample]

test_res <- scores %>%
  group_by(module) %>%
  summarise(
    contrast = "Naive_vs_Healthy",
    delta = mean(score[grp == "Naive"], na.rm = TRUE) - mean(score[grp == "Healthy"], na.rm = TRUE),
    P = tryCatch(wilcox.test(score[grp == "Naive"], score[grp == "Healthy"])$p.value,
                 error = function(e) NA_real_),
    .groups = "drop"
  )
test_res$FDR <- p.adjust(test_res$P, "BH")
print(test_res, n = Inf)
write_csv(test_res, "results/33_blood_exosome/33_blood_module_tests.csv")

# ---------------------------------------------------------------------------
# 4. 外泌体 hub 基因逐基因检验（血液可及标志物候选）
# ---------------------------------------------------------------------------
hub_genes <- c("CD9","CD81","ITGB1","SDC4","ITGAV","YBX1","TSG101","PDCD6IP",
               "LGALS3","LGALS1","CD44","ITGAM","ITGB2","SMPD3")
hub_res <- bind_rows(lapply(intersect(hub_genes, rownames(mat)), function(g) {
  v <- mat[g, ]
  h <- v[names(grp)[grp == "Healthy"]]; nv <- v[names(grp)[grp == "Naive"]]
  tibble(gene = g,
         logFC_naive = mean(nv) - mean(h),
         P_naive = tryCatch(wilcox.test(nv, h)$p.value, error = function(e) NA_real_))
}))
hub_res$FDR <- p.adjust(hub_res$P_naive, "BH")
print(hub_res, n = Inf)
write_csv(hub_res, "results/33_blood_exosome/33_blood_hub_genes.csv")

# ---------------------------------------------------------------------------
# 5. 图
# ---------------------------------------------------------------------------
p1 <- test_res %>%
  filter(!is.na(delta)) %>%
  mutate(module = fct_reorder(module, delta)) %>%
  ggplot(aes(delta, module)) +
  geom_col(fill = "#C9828B") +
  labs(title = "Exosome machinery module scores in peripheral blood (GSE143272, drug-naive 34 vs healthy 51)",
       x = "module score difference vs healthy", y = NULL)

p2 <- hub_res %>%
  mutate(sig = FDR < 0.05, gene = fct_reorder(gene, logFC_naive)) %>%
  ggplot(aes(logFC_naive, gene, fill = sig)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#C9828B", "FALSE" = "#B9C4D4")) +
  labs(title = "Exosome hub genes in blood: drug-naive epilepsy vs healthy",
       x = "log2 expression difference", y = NULL, fill = "FDR<0.05")

pdf("results/33_blood_exosome/33_blood_exosome.pdf", width = 10, height = 5)
print(p1); print(p2)
dev.off()

message("完成。判读：")
message("- 模块阴性 = 与主文 bulk 血液阴性一致，如实报告（需分选单核细胞）")
message("- hub 基因中若有 FDR<0.05 的（尤其 CD9/ITGB1），是患者血样验证的直接靶标")
