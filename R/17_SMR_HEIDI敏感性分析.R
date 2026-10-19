# ============================================================================
# Script      : 17_SMR_HEIDI敏感性分析.R
# Purpose     : SMR + HEIDI sensitivity analysis: LGALS3 cortical eQTL x focal-HS epilepsy GWAS (prepares SMR input + batch files).
# Manuscript  : Methods 2.12 / Results 3.13
# Input       : cortical eQTL data (as in 16); epilepsy GWAS; external SMR binary
# Output      : SMR_analysis/ input files + .bat; interpretation via read_smr()
# Dependencies: data.table, dplyr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 17_SMR_HEIDI.R（最终版，整段一次性运行）
## LGALS3 皮层 eQTL × HS 癫痫 GWAS 的 SMR + HEIDI 敏感性分析
## 依赖：脚本16已跑通（生成皮层 eQTL 数据流程一致）
## 产出：SMR_analysis/ 文件夹（输入文件 + bat）→ 手动双击 bat → read_smr() 解读
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

## ---- 0. 配置 ----
## [EN] Step 0: Configuration
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")
SMR_DIR  <- file.path(BASE_DIR, "SMR_analysis")
dir.create(SMR_DIR, showWarnings = FALSE, recursive = TRUE)

GENE_SYMBOL <- "LGALS3"
GENE_CHR    <- 14
## MetaBrain 皮层 chr14 文件（与脚本16相同）
CORTEX_FILE <- file.path(BASE_DIR, "2021-07-23-cortex-EUR-80PCs-chr14.txt.gz")

## ---- 1. 读取皮层 eQTL，生成 SMR 格式（全 cis SNP，不做 P 过滤）----
cat("---- 读取 MetaBrain 皮层 eQTL ----\n")
d <- fread(CORTEX_FILE, sep = "\t")
d <- d[toupper(d$GeneSymbol) == toupper(GENE_SYMBOL), ]
if (!nrow(d)) stop("皮层文件中未找到 ", GENE_SYMBOL)

f4  <- tstrsplit(d$SNP, ":")          # 14:54673482:rs2057370:C_T
a12 <- tstrsplit(f4[[4]], "_")

eqtl_smr <- data.frame(
  SNP       = f4[[3]],
  Chr       = as.integer(f4[[1]]),
  BP        = as.integer(f4[[2]]),
  A1        = a12[[1]],
  A2        = a12[[2]],
  Freq      = d$SNPEffectAlleleFreq,
  Probe     = GENE_SYMBOL,
  Probe_Chr = GENE_CHR,
  Probe_bp  = as.integer(d$GenePos[1]),
  Gene      = GENE_SYMBOL,
  Orientation = "+",
  b         = d$MetaBeta,
  se        = d$MetaSE,
  p         = d$MetaP,
  stringsAsFactors = FALSE
)
eqtl_smr <- eqtl_smr[complete.cases(eqtl_smr), ]
cat("SMR eQTL 输入 SNP 数:", nrow(eqtl_smr), "\n")
write.table(eqtl_smr, file.path(SMR_DIR, "LGALS3_cortex_eqtl.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 2. 提取 HS GWAS（ieu-b-14）并转 .ma 格式 ----
cat("---- 提取 ieu-b-14 GWAS ----\n")
gwas <- ieugwasr::associations(eqtl_smr$SNP, "ieu-b-14")
ma <- data.frame(
  SNP  = gwas$rsid,
  A1   = gwas$ea,
  A2   = gwas$nea,
  freq = ifelse(is.na(gwas$eaf), 0.5, gwas$eaf),
  b    = gwas$beta,
  se   = gwas$se,
  p    = gwas$p,
  N    = 30480,
  stringsAsFactors = FALSE
)
ma <- ma[complete.cases(ma), ]
cat("GWAS .ma SNP 数:", nrow(ma), "\n")
write.table(ma, file.path(SMR_DIR, "HS_gwas.ma"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 3. 写出 run_SMR.bat（已修正路径与 qfile 转换）----
bat <- c(
  "@echo off",
  "cd /d %~dp0",
  "set SMR=..\\smr-1.3.1-win-x86_64\\smr-1.3.1-win-x86_64\\smr-1.3.1-win.exe",
  "set REF=..\\1000G_EUR\\g1000_eur",
  "",
  "rem ==== 第1步：文本 eQTL 转二进制 BESD（qfile 格式）====",
  "%SMR% --qfile LGALS3_cortex_eqtl.txt --make-besd --out LGALS3_cortex_eqtl",
  "",
  "rem ==== 第2步：SMR + HEIDI ====",
 " %SMR% --bfile %REF% --gwas-summary HS_gwas.ma --beqtl-summary LGALS3_cortex_eqtl ^
    --out LGALS3_HS_SMR_multi --thread-num 4 --diff-freq-prop 0.9 --smr-multi",
  "  --out LGALS3_HS_SMR --thread-num 4 --diff-freq-prop 0.9",
  "",
  "echo.",
  "echo 完成，结果见 LGALS3_HS_SMR.smr",
  "pause"
)
writeLines(bat, file.path(SMR_DIR, "run_SMR.bat"))
cat("bat 已写出:", file.path(SMR_DIR, "run_SMR.bat"), "\n")
cat("\n>>> 请到 SMR_analysis 文件夹双击 run_SMR.bat，跑完回来执行 read_smr() <<<\n\n")

## ---- 4. SMR 结果读取与解读（bat 跑完后执行本行）----
read_smr <- function() {
  f <- file.path(SMR_DIR, "LGALS3_HS_SMR.smr")
  if (!file.exists(f)) stop("未找到 ", f, "，请先双击运行 run_SMR.bat")
  res <- read.delim(f)
  print(res)
  cat("\n========== 关键解读 ==========\n")
  cat("b_SMR   =", round(res$b_SMR[1], 5),
      "（皮层 eQTL MR 方向 b=-0.0158，预期同号为负）\n")
  cat("p_SMR   =", signif(res$p_SMR[1], 3),
      ifelse(res$p_SMR[1] < 0.05, "  → 显著", "  → 不显著"), "\n")
  cat("p_HEIDI =", signif(res$p_HEIDI[1], 3),
      ifelse(res$p_HEIDI[1] > 0.05,
             "  → 通过：信号来自同一因果变异（非LD混杂），可回应coloc欠功效的质疑",
             "  → 未通过：存在LD混杂，结论需谨慎表述"), "\n")
  invisible(res)
}
read_smr <- function() {
  if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
  f <- file.path(BASE_DIR, "SMR_analysis", "LGALS3_HS_SMR.smr")
  if (!file.exists(f)) stop("未找到 ", f, "，请先双击运行 run_SMR.bat")
  res <- read.delim(f)
  print(res)
  cat("\n========== 关键解读 ==========\n")
  cat("b_SMR   =", round(res$b_SMR[1], 5),
      "（皮层 eQTL MR 方向 b=-0.0158，预期同号为负）\n")
  cat("p_SMR   =", signif(res$p_SMR[1], 3),
      ifelse(res$p_SMR[1] < 0.05, "  → 显著", "  → 不显著"), "\n")
  cat("p_HEIDI =", signif(res$p_HEIDI[1], 3),
      ifelse(res$p_HEIDI[1] > 0.05,
             "  → 通过：信号来自同一因果变异（非LD混杂），可回应coloc欠功效的质疑",
             "  → 未通过：存在LD混杂，结论需谨慎表述"), "\n")
  invisible(res)
}
read_smr()
res2 <- read.delim(file.path(BASE_DIR, "SMR_analysis", "LGALS3_HS_SMR_multi.msmr"))
print(res2)
