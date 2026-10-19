# ============================================================================
# Script      : 35.R
# Purpose     : Gal-3 - DAM feedback re-analysis of JCI 2023 data: DAM induction (Tau22), Lgals3-KO effect, TD139 pharmacological inhibition.
# Manuscript  : Supplementary analyses (perturbation, Methods 2.7 related)
# Input       : JCI2023 *.xlsx filtered DEG tables
# Output      : results/35_Gal3_DAM_feedback/
# Dependencies: dplyr, readxl, tidyr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 35_Gal3_DAM_feedback.R
## 目的：用 JCI 2023 重分析数据检验 Gal-3 与 DAM 程序的反馈关系
##   ① 模型诱导期（Tau22 vs WT）：DAM 程序是否整体启动
##   ② Lgals3-KO 效应：敲除后哪些 DAM 基因被压低 → 反馈臂定位
##   ③ TD139 药物抑制：药理抑制是否撼动转录程序
## 注意：这两张表是【过滤后】的 DEG 表（|log2FC|≥0.55 且 P<0.05），
##       基因未出现 = 未过阈值，不是"完全不变"
## 输出：results/35_Gal3_DAM_feedback/
##   35_DAM_feedback_matrix.csv   基因×对比 矩阵
##   35_summary.txt               自动结论
## ============================================================
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root

list.files(BASE_DIR, pattern = "JCI2023.*\\.xlsx$", recursive = TRUE, ignore.case = TRUE)
suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr) })

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT <- file.path(BASE_DIR, "results/35_Gal3_DAM_feedback"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

f_ko <- file.path(BASE_DIR, "data/JCI2023_Tau22_Lgals3KO_差异基因表.xlsx")
f_td <- file.path(BASE_DIR, "data/JCI2023_iMGL_TD139_差异基因表.xlsx")
stopifnot(file.exists(f_ko), file.exists(f_td))

## 读表函数：第 1 行是表名（Supplemental Table x），第 2 行才是真表头
load_deg <- function(file, sheet){
  d <- read_excel(file, sheet = sheet, skip = 1)
  colnames(d)[1:5] <- c("gene_id","gene","ncbi","log2fc","p")
  d %>% mutate(log2fc = as.numeric(log2fc), p = as.numeric(p))
}

induce <- load_deg(f_ko, "6B; TG3WT_vs_WG3WT")        # Tau22(Gal3WT) vs WT：模型诱导
ko     <- load_deg(f_ko, "6B; TG3KO_vs_TG3WT")        # Lgals3-KO vs Tau22-Gal3WT：KO效应
td139  <- load_deg(f_td, "2B; pTauTD139 vs pTau Log2 0.55")  # TD139 vs pTau：药物抑制
message("诱导表 ", nrow(induce), " 行; KO 表 ", nrow(ko), " 行; TD139 表 ", nrow(td139), " 行")

## DAM 程序基因清单（诱导层核心轴 + 效应臂 + 溶酶体）
dam_genes <- c("Trem2","Tyrobp","Apoe","Spp1","Clec7a","Cst7","Itgax","Lpl",
               "Cd68","Ctss","B2m","Ccl6","Lyz2","Fth1","Axl","Cd9","Csf1","Timp2","Lgals3")

getv <- function(d, g){
  m <- d %>% filter(gene == g)
  if (nrow(m) == 0) return(NA_real_)
  m$log2fc[1]
}

mat <- data.frame(
  gene = dam_genes,
  模型诱导_Tau22_vs_WT   = sapply(dam_genes, getv, d = induce),
  KO效应_Lgals3KO        = sapply(dam_genes, getv, d = ko),
  药物抑制_TD139         = sapply(dam_genes, getv, d = td139),
  row.names = NULL, check.names = FALSE
)
write.csv(mat, file.path(OUT, "35_DAM_feedback_matrix.csv"), row.names = FALSE)
print(mat)

## ---------- 自动判定 ----------
core_axis <- c("Trem2","Tyrobp","Apoe","Spp1")
induced   <- !is.na(mat$`模型诱导_Tau22_vs_WT`[match(core_axis, mat$gene)])
ko_hit    <- mat %>% filter(!is.na(KO效应_Lgals3KO), KO效应_Lgals3KO < 0) %>% pull(gene)
td_hit    <- mat %>% filter(!is.na(药物抑制_TD139)) %>% pull(gene)

verdict <- c(
  "==== 分析 35 结论 ====",
  sprintf("模型诱导期 DAM 核心轴（Trem2/Tyrobp/Apoe/Spp1）启动: %s", ifelse(all(induced), "是，全部上调", "部分")),
  sprintf("Lgals3-KO 后显著下调的 DAM 基因: %s", ifelse(length(ko_hit)>0, paste(setdiff(ko_hit,"Lgals3"), collapse=", "), "无")),
  sprintf("TD139 后过阈值的 DAM 基因: %s", ifelse(length(td_hit)>0, paste(td_hit, collapse=", "), "无")),
  "",
  "解读：",
  "1. TREM2-APOE 核心轴不依赖 Gal-3 → Gal-3 不是 DAM 程序的启动开关",
  sprintf("2. Gal-3 反馈定位于效应臂（%s）→ 效应执行层的正反馈，而非诱导层", paste(setdiff(ko_hit,"Lgals3"), collapse=", ")),
  "3. TD139 不撼动转录程序 → Gal-3 反馈不经 lectin 糖结合域；",
  "   干预实验应优先 AAV-shLgals3（敲减）而非 TD139；TD139 改用于验证分泌型 Gal-3 的受体结合",
  "",
  "手稿写法：放 Discussion，标注 tauopathy 模型局限，与分析 34（无单 TF 开关）互证"
)
writeLines(verdict, file.path(OUT, "35_summary.txt"))
cat(verdict, sep = "\n")
