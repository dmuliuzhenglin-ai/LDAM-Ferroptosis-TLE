# ============================================================================
# Script      : 38_LGALS3_coloc_SMR.R
# Purpose     : Colocalization (coloc.abf) and SMR/HEIDI between LGALS3 and epilepsy GWAS to distinguish linkage from shared causal variants.
# Manuscript  : Methods 2.10 / Results 3.11
# Input       : LGALS3 cis-eQTL; epilepsy GWAS summary stats
# Output      : results/38_LGALS3_coloc/
# Dependencies: coloc, data.table, ieugwasr, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 38_LGALS3_coloc_SMR.R
## 目的：LGALS3 与癫痫 GWAS 的共定位（coloc）+ SMR/HEIDI
##   —— 回应"MR 关联是连锁还是因果共享变异"这一审稿人必问点
## 判定标准：
##   coloc  PP.H4 > 0.7 → 共享因果变异（支持 Gal-3 因果）
##          PP.H3 > 0.7 → 两个独立信号（MR 关联是 LD 假象）
##          均不高     → 效能不足，如实报告
##   SMR    p_SMR < 0.05 且 p_HEIDI > 0.01 → 通过
## 输出：results/38_LGALS3_coloc/
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse); library(data.table)
})
## 首次运行先装：install.packages("coloc")
suppressPackageStartupMessages(library(coloc))

set.seed(2024)
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/38_LGALS3_coloc", showWarnings = FALSE, recursive = TRUE)

## ================= 参数区（按需修改） =================
## [EN] Parameter block (edit as needed)
## LGALS3 位置（GRCh37/hg19，与 OpenGWAS 多数数据一致）
CHR     <- 14
GENE_LO <- 55612775
GENE_HI <- 55629907
WINDOW  <- 5e5                       # cis 窗口 ±500 kb
WIN_LO  <- GENE_LO - WINDOW
WIN_HI  <- GENE_HI + WINDOW

## OpenGWAS ID（如自动探测失败，手工填入你 MR 时用过的 ID）
EQTL_ID <- "eqtl-a-ENSG00000131981"  # eQTLGen 血液 cis-eQTL, LGALS3
GWAS_ID <- ""                        # ← 留空则自动从旧 MR 文件探测；建议填海马硬化型 focal epilepsy 的 ID

## ================= 1. 数据获取 =================
## [EN] Step 1: data acquisition
## 策略 A：复用本地 MR 时的 harmonised 数据（最省事、口径完全一致）
find_local <- function(){
  fs <- list.files("results", pattern = "\\.(csv|tsv|rds|txt)$",
                   recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  fs[grepl("harmonis|harmoniz|coloc|eqtl.*lgals3|lgals3.*eqtl|MR_", basename(fs), ignore.case = TRUE)]
}
local_files <- find_local()
message("本地候选文件: ", length(local_files))
if (length(local_files)) print(basename(local_files))

dat <- NULL
if (length(local_files) >= 1) {
  f <- local_files[1]
  message("尝试读取: ", f)
  dat <- tryCatch({
    if (grepl("\\.rds$", f, ignore.case = TRUE)) readRDS(f) else fread(f) |> as.data.frame()
  }, error = function(e) NULL)
}

## 策略 B：OpenGWAS 在线拉取（需要 .Renviron 里有效 JWT）
if (is.null(dat)) {
  message("本地无可用文件 → 走 OpenGWAS 在线模式")
  suppressPackageStartupMessages(library(ieugwasr))
  if (GWAS_ID == "") stop("请把 MR 时用的 GWAS ID 填到脚本顶部 GWAS_ID 参数（海马硬化型 focal epilepsy 那个）")

  eqtl <- ieugwasr::associations(variants = NULL, id = EQTL_ID)  # 占位，实际用下述区域提取
  ## 区域提取（ieugwasr 无区域接口时用 API 直接拉）
  tok <- Sys.getenv("OPENGWAS_JWT")
  if (tok == "") stop("OpenGWAS JWT 为空：请到 https://api.opengwas.io 后台重新生成，写入 .Renviron 的 OPENGWAS_JWT，然后重启 R")

  pull_region <- function(id){
    url <- sprintf("https://api.opengwas.io/api/associations/%s?chr=%d&start=%d&stop=%d",
                   id, CHR, WIN_LO, WIN_HI)
    resp <- httr::GET(url, httr::add_headers(Authorization = paste("Bearer", tok)))
    httr::stop_for_status(resp)
    as.data.frame(httr::content(resp, as = "parsed"))
  }
  e <- pull_region(EQTL_ID)   # trait 1: eQTL
  g <- pull_region(GWAS_ID)   # trait 2: GWAS

  ## 统一列名（OpenGWAS 返回: rsid, chr, position, beta, se, p, n, ea, nea, eaf...）
  std <- function(d, n_col = "n"){
    d$position <- as.numeric(d$position)
    data.frame(rsid = d$rsid, pos = d$position,
               beta = as.numeric(d$beta), se = as.numeric(d$se),
               p = as.numeric(d$p), N = as.numeric(d[[n_col]]),
               stringsAsFactors = FALSE)
  }
  e2 <- std(e); g2 <- std(g)
  dat <- list(eqtl = e2, gwas = g2)
}

## ================= 2. 整理为 coloc 输入 =================
## [EN] Step 2: format coloc input
if (!is.list(dat)) {
  ## 本地 harmonised 单列格式：自动识别 beta/se/p 列
  cn <- names(dat)
  pick <- function(pat) cn[grepl(pat, cn, ignore.case = TRUE)][1]
  mk <- function(prefix){
    data.frame(
      rsid = dat[[pick("^SNP$|rsid")]],
      beta = as.numeric(dat[[pick(paste0("beta\\.", prefix, "$|beta_", prefix))]]),
      se   = as.numeric(dat[[pick(paste0("se\\.", prefix, "$|se_", prefix))]]),
      p    = as.numeric(dat[[pick(paste0("pval\\.", prefix, "$|p_", prefix))]]),
      stringsAsFactors = FALSE)
  }
  e2 <- mk("exposure"); g2 <- mk("outcome")
  ## 样本量（ harmonised 里常没有，手工设）
  e2$N <- 31684      # eQTLGen 2019 cis-eQTL 样本量
  g2$N <- NA_real_   # 填你 GWAS 的样本量；case/control 比例另填下面 s
  dat <- list(eqtl = e2, gwas = g2)
}

e2 <- dat$eqtl; g2 <- dat$gwas
m <- inner_join(e2, g2, by = "rsid", suffix = c(".e", ".g")) |>
  filter(!is.na(beta.e), !is.na(se.e), !is.na(beta.g), !is.na(se.g))
message("窗口内共同 SNP: ", nrow(m), "（", WIN_LO, "-", WIN_HI, "）")
stopifnot(nrow(m) >= 20)

## ================= 3. coloc.abf =================
D1 <- list(beta = m$beta.e, varbeta = m$se.e^2, snp = m$rsid,
           N = median(m$N.e, na.rm = TRUE), type = "quant")
N_gwas <- suppressWarnings(median(m$N.g, na.rm = TRUE))
if (is.na(N_gwas)) N_gwas <- 50000   # 若 GWAS 无 n 列，按实际填
D2 <- list(beta = m$beta.g, varbeta = m$se.g^2, snp = m$rsid,
           N = N_gwas, type = "cc", s = 0.3)   # s = 病例比例，按实际 GWAS 改

cc <- coloc.abf(D1, D2)
pp <- cc$summary
print(pp)

## 后验概率解读
h4 <- pp["PP.H4.abf"]; h3 <- pp["PP.H3.abf"]
verdict <- if (h4 > 0.7) {
  sprintf("【共享因果变异】PP.H4 = %.3f：LGALS3 eQTL 与癫痫 GWAS 信号共定位，MR 关联支持因果解释", h4)
} else if (h3 > 0.7) {
  sprintf("【独立信号】PP.H3 = %.3f：两个信号是 LD 假象，MR 结论需降级为'遗传关联'", h3)
} else {
  sprintf("【效能不足】PP.H3 = %.3f / PP.H4 = %.3f：无法区分，建议如实报告并注明 coloc 不结论", h3, h4)
}

## 区域 top SNP（供区域堆叠图/LocusZoom 描述）
top <- m |> arrange(p.g) |> slice_head(n = 5) |>
  dplyr::select(rsid, beta.e, p.e, beta.g, p.g)

## ================= 4. SMR / HEIDI（可选，需要 SMR 二进制） =================
## [EN] Step 4: SMR / HEIDI (optional; needs SMR binary)
## SMR 是独立命令行工具（Yang lab），R 里无法直接跑：
##   1) 下载 smr（https://yanglab.westlake.edu.cn/software/smr）
##   2) 准备：eQTL 的 .besd 文件 + GWAS sumstats（SMR 格式）+ 1000G EUR LD 参考（.bed/.bim/.fam）
##   3) 命令：
##      smr --bfile 1000G_EUR --gwas-summary gwas.ma --beqtl-summary lgals3_eqtl \
##          --out smr_lgals3 --thread-num 4 --diff-freq-prop 0.5
##   判定：p_SMR < 0.05 且 p_HEIDI > 0.01 才通过
## 如果你之前 MR 用的是 MetaBrain 皮层 eQTL，SMR 用同一份 eQTL 口径最一致。

## ================= 5. 写出 =================
## [EN] Step 5: Write outputs
out <- "results/38_LGALS3_coloc"
write.csv(m, file.path(out, "38_window_snps.csv"), row.names = FALSE)
write.csv(as.data.frame(t(pp)), file.path(out, "38_coloc_posteriors.csv"), row.names = FALSE)
txt <- c("==== 分析 38：LGALS3 coloc ====",
         sprintf("cis 窗口 chr%d:%d-%d，共同 SNP %d 个", CHR, WIN_LO, WIN_HI, nrow(m)),
         sprintf("PP.H0=%.3g PP.H1=%.3g PP.H2=%.3g PP.H3=%.3g PP.H4=%.3g",
                 pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"], h3, h4),
         "", verdict, "",
         "GWAS 信号 top 5 SNP:", capture.output(print(top, row.names = FALSE)))
cat(txt, file = file.path(out, "38_summary.txt"), sep = "\n")
cat(txt, sep = "\n")
