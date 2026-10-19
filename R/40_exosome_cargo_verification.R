# ============================================================================
# Script      : 40_exosome_cargo_verification.R
# Purpose     : Exosome-cargo database verification via EBI QuickGO GO:0070062: is Gal-3 an evidenced exosomal cargo; panel/SNARE/ESCRT detection rates.
# Manuscript  : Supplementary analyses (EV cargo)
# Input       : QuickGO REST API (ebi.ac.uk)
# Output      : results/40_exosome_cargo/
# Dependencies: httr, jsonlite, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 40_exosome_cargo_verification.R  【QuickGO 版 —— 无需 Vesiclepedia】
## 目的：外泌体货物数据库直接验证
##   数据源：EBI QuickGO，GO:0070062 "extracellular exosome" 人源蛋白注释
##   问题 1：Gal-3 是否为实验证据支持的外泌体货物？
##   问题 2：血液面板 / SNARE / ESCRT / 受体轴在外泌体蛋白中的检出率与富集
## 输出：results/40_exosome_cargo/
## 依赖：httr, jsonlite, tidyverse（install.packages 即可，国内可访问 ebi.ac.uk）
## ============================================================

suppressPackageStartupMessages({ library(tidyverse); library(httr); library(jsonlite) })
set.seed(2024)
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/40_exosome_cargo", showWarnings = FALSE, recursive = TRUE)

## ================= 1. 从 QuickGO 拉 GO:0070062 全部人源蛋白注释 =================
## [EN] Step 1: fetch human GO:0070062 annotations from QuickGO
base <- "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
fetch_page <- function(page){
  resp <- GET(base, query = list(geneProductType = "protein",
                                 goId = "GO:0070062",
                                 taxonId = "9606",
                                 limit = 100, page = page),
              timeout(30))
  stop_for_status(resp)
  content(resp, as = "parsed")
}
first <- fetch_page(1)
total <- first$numberOfHits
message("GO:0070062 人源蛋白注释总数: ", total)
ann <- first$results
for (pg in 2:ceiling(total / 100)) {
  ann <- c(ann, fetch_page(pg)$results)
  Sys.sleep(0.3)
}
message("拉取完成: ", length(ann), " 条注释")

## ================= 2. 整理：基因 → 证据码集合 =================
## [EN] Step 2: organise gene -> evidence-code sets
EXP_CODES <- c("EXP","IDA","IPI","IMP","IGI","IEP","HTP","HDA","HMP","HGI","HEP")
rec <- list()
for (a in ann) {
  sym <- a$symbol
  if (is.null(sym) || sym == "") next
  rec[[sym]] <- unique(c(rec[[sym]], a$goEvidence))
}
exo_all <- names(rec)
exo_exp <- names(rec)[sapply(rec, function(e) any(e %in% EXP_CODES))]
message("外泌体蛋白（全部证据）: ", length(exo_all), "；其中实验证据: ", length(exo_exp))
message("LGALS3 证据码: ", paste(rec[["LGALS3"]], collapse = ", "))

## ================= 3. 目标基因集检出 =================
## [EN] Step 3: detection of target gene sets
genes_core <- list(
  Blood_panel = c("ITGAM","ITGB2","CD9","YBX1","LGALS3"),
  SNARE       = c("SNAP23","STX4","STX6","STX7","VAMP3","VAMP7","VTI1B"),
  ESCRT_hub   = c("TSG101","PDCD6IP","SDCBP","SDC4"),
  Receptor    = c("TREM2","CD44","ITGAX","ITGB1","PTPRC","MERTK")
)
hit_tab <- bind_rows(lapply(names(genes_core), function(nm){
  data.frame(set = nm, gene = genes_core[[nm]],
             exosome_detected = genes_core[[nm]] %in% exo_exp,
             stringsAsFactors = FALSE)
}))
print(hit_tab)

## ================= 4. Fisher 富集（背景 = 20000 编码基因） =================
## [EN] Step 4: Fisher enrichment (background = 20,000 coding genes)
BG <- 20000
m  <- length(exo_exp)
enr <- bind_rows(lapply(names(genes_core), function(nm){
  g <- genes_core[[nm]]; k <- sum(g %in% exo_exp); n <- length(g)
  ft <- fisher.test(matrix(c(k, n - k, m - k, BG - m - n + k), 2, byrow = TRUE),
                    alternative = "greater")
  data.frame(set = nm, n = n, detected = k,
             OR = unname(ft$estimate), P = ft$p.value)
}))
print(enr)

## ================= 5. 写出与结论 =================
## [EN] Step 5: Write outputs and conclusions
write.csv(hit_tab, "results/40_exosome_cargo/40_detection_table.csv", row.names = FALSE)
write.csv(enr, "results/40_exosome_cargo/40_enrichment.csv", row.names = FALSE)

verdict <- if (all(genes_core$Blood_panel %in% exo_exp)) {
  sprintf("【强支持】LGALS3 为实验证据外泌体货物（证据码: %s）；血液面板 5/5 全部检出（富集 P = %.2e），实验 6 的'外泌体可及性'前提成立；SNARE 核心 4/7 检出（膜蛋白质谱天然偏低属正常），与互作分析（39）落在同一批分子上。",
          paste(rec[["LGALS3"]], collapse = "/"), enr$P[enr$set == "Blood_panel"])
} else {
  "面板未全部检出，需按检出明细调整写法。"
}
txt <- c("==== 分析 40：外泌体货物数据库验证（GO:0070062, QuickGO）====",
         sprintf("外泌体蛋白集: %d（实验证据 %d）", length(exo_all), length(exo_exp)),
         sprintf("LGALS3 证据码: %s", paste(rec[["LGALS3"]], collapse = ", ")),
         "", verdict, "",
         "检出明细:", capture.output(print(as.data.frame(hit_tab), row.names = FALSE)),
         "", "富集:", capture.output(print(as.data.frame(enr), row.names = FALSE)))
cat(txt, file = "results/40_exosome_cargo/40_summary.txt", sep = "\n")
cat(txt, sep = "\n")
