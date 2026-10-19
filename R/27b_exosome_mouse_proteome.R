#!/usr/bin/env Rscript
# ============================================================================
# Script      : 27b_exosome_mouse_proteome.R
# Purpose     : External independent validation: chronic mouse TLE hippocampal proteome (KA 35d LFQ, 3v3) exosome-module Fisher enrichment.
# Manuscript  : Supplementary analyses (cross-species proteome)
# Input       : 12974_2023_2840_MOESM3_ESM.xlsx (published supplementary table)
# Output      : results/27b_exosome_mouse_proteome/ 4 CSV
# Dependencies: readxl, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================
# ===========================================================================
# 27b_exosome_mouse_proteome.R
# 第四部分·外部独立验证：慢性期小鼠 TLE 海马蛋白组（KA 35 天，3v3，LFQ）
# 来源：文献补充表 MOESM3（J Neuroinflammation 2023, 文章号 2840，Table S2）
#       课题组引用该文，使用时正文须规范引用
# 内容：
#   (1) 读入 6085 蛋白全量表（含 LFQ 强度、TLE/Con Ratio、P value）
#   (2) 外泌体机器 7 模块 Fisher 富集（整体 + 上调/下调方向分层）
#   (3) 逐基因关键报告（LGALS3 预期全 0 = 未检出）
# 输出：results/27b_exosome_mouse_proteome/ 下 4 个 CSV
# 依赖：tidyverse, readxl
# ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/27b_exosome_mouse_proteome", showWarnings = FALSE, recursive = TRUE)

find1 <- function(pattern, exclude = NULL) {
  f <- list.files("data", pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!is.null(exclude)) f <- f[!grepl(exclude, basename(f))]
  if (length(f) == 0) stop("找不到匹配文件: ", pattern, "（把 MOESM3 xlsx 放进 data/ 文件夹）")
  message("使用文件: ", f[1])
  f[1]
}

# ---------------------------------------------------------------------------
# 1. 读表：全量蛋白表在 sheet "Protein  (2)"（6085 行）
#    若该 sheet 名在你本地文件里不同，先用 excel_sheets(f) 打印核对
# ---------------------------------------------------------------------------
f <- find1("12974_2023_2840.*MOESM3.*xlsx")
print(excel_sheets(f))
raw <- read_excel(f, sheet = "Protein  (2)")   # 注意 sheet 名以打印结果为准
stopifnot(all(c("Gene name", "TLE/Con Ratio", "TLE/Con P value") %in% colnames(raw)))

df <- raw %>%
  filter(!is.na(`Gene name`)) %>%
  mutate(gene = toupper(str_remove(`Gene name`, ";.*"))) %>%
  mutate(sig    = !is.na(`TLE/Con P value`) & `TLE/Con P value` < 0.05,
         sig_up = sig & !is.na(`TLE/Con Ratio`) & `TLE/Con Ratio` > 1,
         sig_dn = sig & !is.na(`TLE/Con Ratio`) & `TLE/Con Ratio` < 1)

message("检出蛋白: ", nrow(df), " | 显著(P<0.05): ", sum(df$sig))
message("LGALS3 检出: ", ifelse("LGALS3" %in% df$gene,
        sum(df[df$gene=="LGALS3", grepl("LFQ", colnames(df))] > 0, na.rm = TRUE), "未检出"))

# ---------------------------------------------------------------------------
# 2. 外泌体机器 7 模块 Fisher 富集（universe = 全部 6085 检出蛋白）
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

fisher_one <- function(gset, cond) {
  a <- sum(df$gene %in% gset & cond)          # 模块内&条件
  b <- sum(df$gene %in% gset & !cond)
  c <- sum(!(df$gene %in% gset) & cond)
  d <- sum(!(df$gene %in% gset) & !cond)
  ft <- fisher.test(matrix(c(a, b, c, d), 2, 2))
  tibble(n_sig_dir = a, OR = unname(ft$estimate), P = ft$p.value)
}

res_overall <- bind_rows(lapply(names(gene_sets), function(m)
  fisher_one(gene_sets[[m]], df$sig) %>% mutate(module = m, direction = "any"))) %>%
  select(module, direction, everything())
res_overall$FDR <- p.adjust(res_overall$P, "BH")

res_dir <- bind_rows(lapply(names(gene_sets), function(m) bind_rows(
  fisher_one(gene_sets[[m]], df$sig_up) %>% mutate(module = m, direction = "up"),
  fisher_one(gene_sets[[m]], df$sig_dn) %>% mutate(module = m, direction = "down")))) %>%
  select(module, direction, everything())
res_dir$FDR <- p.adjust(res_dir$P, "BH")

# ---------------------------------------------------------------------------
# 3. 关键基因逐基因报告
# ---------------------------------------------------------------------------
key_genes <- c("CD9","CD81","ICAM1","ITGB1","ITGAV","SDC1","SDC4","HSPG2",
               "YBX1","SYNCRIP","HNRNPA2B1","FMR1",
               "TSG101","PDCD6IP","CHMP4B",
               "RAB27A","RAB27B","SMPD3","LGALS3","GFAP","GPX4","ACSL4","PLIN2")
key_tbl <- df %>% filter(gene %in% key_genes) %>%
  select(gene, ratio = `TLE/Con Ratio`, P = `TLE/Con P value`) %>%
  arrange(desc(ratio))

# ---------------------------------------------------------------------------
# 4. 输出
# ---------------------------------------------------------------------------
write_csv(res_overall,  "results/27b_exosome_mouse_proteome/27b_module_enrichment_overall.csv")
write_csv(res_dir,      "results/27b_exosome_mouse_proteome/27b_module_enrichment_directional.csv")
write_csv(key_tbl,      "results/27b_exosome_mouse_proteome/27b_key_genes.csv")
write_csv(df %>% select(gene, ratio = `TLE/Con Ratio`, P = `TLE/Con P value`),
          "results/27b_exosome_mouse_proteome/27b_full_with_stats.csv")

print(res_dir %>% arrange(P))
print(key_tbl, n = Inf)

message("完成。判读要点：")
message("- 预期阳性结果：Exosome_uptake 模块 'up' 方向 FDR<0.05（CD9/ITGB1/SDC4/ITGAV/CD81 全线上调）")
message("- 预期阴性：Rab_secretory 不动、SMPD3 下调、LGALS3 未检出 —— 均须如实报告")
message("- 该队列定位：跨物种外泌体机器激活的独立证据；'装载驱动 Gal-3' 仍由大鼠自有蛋白组支撑")
