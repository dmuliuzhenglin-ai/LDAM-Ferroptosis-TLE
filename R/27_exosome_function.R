#!/usr/bin/env Rscript
# ============================================================================
# Script      : 27_exosome_function.R
# Purpose     : Exosome-machinery (uptake/loading/secretion) analysis across human chronic TLE (GSE168375), mouse KA time course (GSE88992) and GSE190451.
# Manuscript  : Supplementary analyses (extracellular-vesicle modules)
# Input       : GSE168375 / GSE88992 / GSE190451 matrices; built-in module gene sets
# Output      : results/27_* tables/figures
# Dependencies: ggplot2, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================
# ===========================================================================
# 27_exosome_function.R  【v2 完整修正版】
# 第四部分：外泌体功能生信分析
# 问题：LGALS3 蛋白在慢性期持续升高（5.1x）而 mRNA 不动（0.605, P=0.163），
#       是否由外泌体装载/分泌环节驱动？
# 队列：GSE168375（人慢性期 31v12）+ GSE88992（小鼠 KA 6/12/24h 时序）
#       + GSE190451（人验证 3v3）
# v2 修复：
#   - to_mouse：首字母大写（TSG101 -> Tsg101）
#   - GSE88992：时间点配对对照（KA_6h vs Con_6h），对照池不再虚构 Con_0
#   - GSE190451：分组硬编码（sample/group 两列），stopifnot 保险丝
#   - pivot_wider 前先 select 三列，防 n_genes 串行
#   - 相关分析加 100 组随机基因集置换背景（否则全局共表达伪显著）
#   - 画图前 filter(!is.na(...))
# 依赖：tidyverse, ggplot2
# ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/27_exosome_function", showWarnings = FALSE, recursive = TRUE)

find1 <- function(pattern, exclude = NULL) {
  f <- list.files("data", pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!is.null(exclude)) f <- f[!grepl(exclude, basename(f))]
  if (length(f) == 0) stop("找不到匹配文件: ", pattern)
  message("使用文件: ", f[1])
  f[1]
}

# ---------------------------------------------------------------------------
# 1. 外泌体功能基因集（人类符号；小鼠矩阵用首字母大写规则匹配）
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
  Exosome_uptake    = c("CD9","CD81","ICAM1","ITGB1","ITGAV","SDC1","SDC4","HSPG2")
)

ref_sets <- list(
  LDAM        = c("PLIN2","DGAT1","DGAT2","ACAT1","ACSL1","CD36","FABP3","HILPDA","SOAT1","LPIN1"),
  Ferroptosis = c("GPX4","ACSL4","LPCAT3","ALOX15","AIFM2","TFRC","FTH1","FTL","SLC7A11","NCOA4")
)

# 人类符号 -> 小鼠符号：首字母大写、其余小写（TSG101 -> Tsg101）
to_mouse <- function(s) paste0(toupper(substr(s, 1, 1)), tolower(substr(s, 2, nchar(s))))

# ---------------------------------------------------------------------------
# 2. 通用函数
# ---------------------------------------------------------------------------
module_scores <- function(expr, symbols_col, sets, mouse = FALSE) {
  mat <- as.matrix(expr[, setdiff(colnames(expr), symbols_col)])
  rownames(mat) <- expr[[symbols_col]]
  if (mouse) sets <- lapply(sets, to_mouse)
  res <- lapply(names(sets), function(nm) {
    g <- intersect(sets[[nm]], rownames(mat))
    if (length(g) < 3)
      return(tibble(sample = colnames(mat), module = nm, score = NA_real_, n_genes = length(g)))
    z <- t(scale(t(mat[g, , drop = FALSE])))
    z[is.na(z)] <- 0
    tibble(sample = colnames(mat), module = nm, score = colMeans(z), n_genes = length(g))
  })
  bind_rows(res)
}

wilcox_module <- function(scores_df, group) {
  scores_df %>% group_by(module) %>% group_modify(function(d, ...) {
    x <- d$score[group[d$sample] == "Case"]; y <- d$score[group[d$sample] == "Con"]
    if (length(x) < 3 || length(y) < 3 || all(is.na(x)))
      return(tibble(logFC = NA_real_, P = NA_real_))
    tibble(logFC = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
           P = tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_))
  }) %>% ungroup()
}

# ---------------------------------------------------------------------------
# 3. GSE168375（人，慢性期，星形胶质细胞富集皮层）——主队列
# ---------------------------------------------------------------------------
message("===== GSE168375 =====")
f168 <- find1("GSE168375_expr\\.csv$")
e168 <- read_csv(f168, show_col_types = FALSE)

# 样本列顺序即 31 TLE + 12 Con（脚本 24 已 stopifnot 验证过该矩阵）
smp168 <- setdiff(colnames(e168), "gene")
stopifnot(length(smp168) == 43)
grp168 <- set_names(c(rep("Case", 31), rep("Con", 12)), smp168)

scores168 <- module_scores(e168, "gene", c(gene_sets, ref_sets), mouse = FALSE)
test168 <- wilcox_module(scores168, grp168) %>% mutate(dataset = "GSE168375_human_chronic")
print(test168, n = Inf)

# ---------------------------------------------------------------------------
# 3b. 样本级相关 + 随机基因集置换背景（100 次）
# ---------------------------------------------------------------------------
message("===== 相关分析（含置换背景）=====")
set.seed(42)
wide168 <- scores168 %>%
  select(sample, module, score) %>%
  pivot_wider(names_from = module, values_from = score)

cor168 <- bind_rows(lapply(names(ref_sets), function(ref) {
  tibble(dataset = "GSE168375", pair = paste0("Exosome_modules_vs_", ref),
         module = names(gene_sets),
         cor = map_dbl(names(gene_sets), ~ tryCatch(
           cor(wide168[[.x]], wide168[[ref]], use = "complete.obs"),
           error = function(e) NA_real_)),
         P = map_dbl(names(gene_sets), ~ tryCatch(
           cor.test(wide168[[.x]], wide168[[ref]])$p.value,
           error = function(e) NA_real_)))
}))

# 置换背景：从表达矩阵随机抽等长基因集，算与 LDAM/Ferroptosis 评分的相关分布
mat168 <- as.matrix(e168[, smp168]); rownames(mat168) <- e168$gene
all_genes168 <- rownames(mat168)
score_one_set <- function(g) {
  g <- intersect(g, all_genes168)
  z <- t(scale(t(mat168[g, , drop = FALSE]))); z[is.na(z)] <- 0
  colMeans(z)
}
ldam_v <- score_one_set(ref_sets$LDAM)
ferro_v <- score_one_set(ref_sets$Ferroptosis)
set.seed(123)
perm_null <- bind_rows(lapply(1:100, function(i) {
  bind_rows(lapply(gene_sets, function(gs) {
    g_rand <- sample(all_genes168, length(gs))
    tibble(pair = c("vs_LDAM", "vs_Ferroptosis"),
           cor = c(cor(score_one_set(g_rand), ldam_v),
                   cor(score_one_set(g_rand), ferro_v)))
  }))
}))
perm_summary <- perm_null %>% group_by(pair) %>%
  summarise(null_mean = mean(cor), null_sd = sd(cor), null_p95 = quantile(cor, 0.95))
write_csv(perm_null, "results/27_exosome_function/27_correlation_permutation_null.csv")
write_csv(perm_summary, "results/27_exosome_function/27_correlation_permutation_summary.csv")
message("置换背景（100 组随机基因集）:")
print(perm_summary)
message("判读：外泌体模块的 r 须显著高于 null_p95 才算特异性共表达")

# ---------------------------------------------------------------------------
# 4. GSE88992（小鼠 KA 急性时序，芯片已 log2，勿再转换）
#    列名：saline_6h_rep1 / kainate_6h_rep1 ...（24h 的 KA 仅 2 重复）
# ---------------------------------------------------------------------------
message("===== GSE88992 =====")
f889 <- find1("GSE88992_expr_by_gene\\.csv$")
g92 <- read_csv(f889, show_col_types = FALSE)
smp92 <- setdiff(colnames(g92), "sym")
message("样本列: ", paste(smp92, collapse = " | "))

cls <- ifelse(grepl("kainate", smp92, ignore.case = TRUE), "KA", "Con")
tim <- str_extract(smp92, "(6|12|24)h")
grp92 <- paste(cls, tim, sep = "_")
names(grp92) <- smp92
print(table(grp92))
stopifnot(sum(grp92 == "Con_6h") == 3, sum(grp92 == "Con_12h") == 3,
          sum(grp92 == "Con_24h") == 3, sum(grp92 == "KA_6h") == 3,
          sum(grp92 == "KA_12h") == 3, sum(grp92 == "KA_24h") == 2)

scores92 <- module_scores(g92, "sym", c(gene_sets, ref_sets), mouse = FALSE)
scores92$grp <- grp92[scores92$sample]

# 保险丝：小鼠基因匹配成功检查（n_genes 应 >= 3）
ng92 <- scores92 %>% group_by(module) %>% summarise(n = first(n_genes))
print(ng92, n = Inf)
stopifnot(all(ng92$n >= 3))

# 时间点配对对照：KA_6h vs Con_6h，依此类推
test92 <- bind_rows(lapply(c("6h", "12h", "24h"), function(tp) {
  ka <- scores92 %>% filter(grp == paste0("KA_", tp))
  co <- scores92 %>% filter(grp == paste0("Con_", tp))
  bind_rows(lapply(unique(ka$module), function(m) {
    x <- ka$score[ka$module == m]; y <- co$score[co$module == m]
    tibble(module = m, grp = paste0("KA_", tp),
           logFC = if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
           P = tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_))
  })) %>% mutate(dataset = "GSE88992_mouse_KA")
}))
print(test92, n = Inf)

# ---------------------------------------------------------------------------
# 5. GSE190451（人验证队列，TLE 3 vs Control 3）
# ---------------------------------------------------------------------------
message("===== GSE190451 =====")
e190 <- read_csv("results/02_GSE190451_normalized.csv", show_col_types = FALSE)
m190 <- read_csv("results/02_GSE190451_sample_meta.csv", show_col_types = FALSE)

colnames(e190) <- trimws(colnames(e190))
grp190 <- set_names(ifelse(m190$group == "TLE", "Case", "Con"), m190$sample)
print(table(grp190))
stopifnot(sum(grp190 == "Case") == 3, sum(grp190 == "Con") == 3)

sym190 <- names(e190)[1]
stopifnot(sym190 == "gene")

scores190 <- module_scores(e190, sym190, c(gene_sets, ref_sets), mouse = FALSE)
scores190 <- scores190 %>% filter(sample %in% names(grp190))
stopifnot(length(unique(scores190$sample)) == 6)

test190 <- wilcox_module(scores190, grp190) %>% mutate(dataset = "GSE190451_validation")
print(test190, n = Inf)

# ---------------------------------------------------------------------------
# 6. 汇总输出
# ---------------------------------------------------------------------------
all_tests <- bind_rows(test168, test92, test190)
write_csv(all_tests, "results/27_exosome_function/27_exosome_module_tests.csv")
write_csv(scores168, "results/27_exosome_function/27_exosome_module_scores_GSE168375.csv")
write_csv(cor168, "results/27_exosome_function/27_exosome_score_correlations.csv")

# ---------------------------------------------------------------------------
# 7. PDF：四个面板
# ---------------------------------------------------------------------------
p1 <- test168 %>%
  filter(!is.na(logFC)) %>%
  mutate(sig = ifelse(is.na(P), "", ifelse(P < 0.05, "*", "")),
         module = fct_reorder(module, logFC)) %>%
  ggplot(aes(logFC, module)) +
  geom_col(aes(fill = P < 0.05)) +
  geom_text(aes(label = sig), hjust = -0.2, size = 5) +
  scale_fill_manual(values = c("TRUE" = "#C9828B", "FALSE" = "#B9C4D4")) +
  labs(title = "Exosome machinery modules: human chronic epilepsy (GSE168375, 31v12)",
       x = "mean z-score difference (TLE - control)", y = NULL, fill = "P<0.05")

p2 <- test92 %>%
  filter(!is.na(logFC)) %>%
  mutate(module = fct_reorder(module, logFC, .fun = max)) %>%
  ggplot(aes(logFC, module, color = grp)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("KA_6h" = "#7A9CC4", "KA_12h" = "#0E3A6B", "KA_24h" = "#C9828B")) +
  labs(title = "Acute KA time course (GSE88992): secretion machinery response",
       x = "module score vs time-matched saline control", y = NULL, color = NULL)

p3 <- cor168 %>%
  filter(!is.na(cor)) %>%
  mutate(sig = ifelse(P < 0.05, "*", ""), module = fct_reorder(module, cor)) %>%
  ggplot(aes(cor, module)) +
  geom_col(aes(fill = pair)) +
  geom_text(aes(label = sig), hjust = -0.2, size = 5) +
  facet_wrap(~ pair, nrow = 1) +
  scale_fill_manual(values = c("#0E3A6B", "#D9A441")) +
  labs(title = "Sample-level correlation: exosome modules vs LDAM / Ferroptosis scores (GSE168375)",
       subtitle = "Compare against permutation null (see 27_correlation_permutation_summary.csv)",
       x = "Pearson r", y = NULL) +
  theme(legend.position = "none")

p4 <- test190 %>%
  filter(!is.na(logFC)) %>%
  mutate(module = fct_reorder(module, logFC)) %>%
  ggplot(aes(logFC, module)) +
  geom_col(fill = "#0E3A6B") +
  labs(title = "Validation cohort GSE190451 (TLE 3v3; Wilcoxon floor P = 0.1)",
       x = "mean z-score difference (TLE - control)", y = NULL)

pdf("results/27_exosome_function/27_exosome_function.pdf", width = 11, height = 9)
print(p1); print(p2); print(p3); print(p4)
dev.off()

message("完成。输出: results/27_exosome_function/ 下 6 个文件")
message("判读矩阵：")
message("- GSE168375 模块差异阴性为预期（星形胶质细胞皮层，LDAM/铁死亡参照集亦阴性）")
message("- GSE190451 3v3 的 P 地板是 0.1，看方向一致性（几/9 同向）而非 P 本身")
message("- GSE88992 看 Rab/Cargo_loading 的时序激活窗口")
