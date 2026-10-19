# ============================================================================
# Script      : 16_脑组织eQTL_MR.R
# Purpose     : Brain-tissue eQTL MR for LGALS3: MetaBrain hippocampus/cortex (EUR) exposures vs epilepsy outcomes.
# Manuscript  : Methods 2.10 / Results 3.11
# Input       : MetaBrain eQTL summary stats (local files); OpenGWAS outcomes
# Output      : results/16_* MR tables
# Dependencies: TwoSampleMR, data.table, dplyr, ieugwasr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================

## ============================================================
## 16_脑组织eQTL_MR.R（最终版，整段一次性运行）
## LGALS3 脑组织 eQTL → 癫痫结局 MR
## 暴露：MetaBrain 海马亚组 + 皮层（EUR）
## 结局：ieu-b-8（全部癫痫）/ ieu-b-10（局灶）/ ieu-b-14（局灶伴HS）
## 注意：运行前确认 OPENGWAS_JWT 已在 .Renviron 中（脚本15已配置）
## ============================================================
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# usethis::edit_r_environ()
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(ieugwasr)
})

## ---- 0. 配置 ----
## [EN] Step 0: Configuration
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")

if (Sys.getenv("OPENGWAS_JWT") == "") {
  stop("未检测到 OPENGWAS_JWT，请先在 .Renviron 中配置令牌并重启 R")
}

GENE_SYMBOL <- "LGALS3"

FILES <- list(
  MetaBrain_Hippocampus = file.path(BASE_DIR, "2021-07-23-hippocampus-EUR-30PCs-chr14.txt.gz"),
  MetaBrain_Cortex      = file.path(BASE_DIR, "2021-07-23-cortex-EUR-80PCs-chr14.txt.gz")
)

OUTCOMES <- c(EPI_ALL = "ieu-b-8", EPI_FOCAL = "ieu-b-10", EPI_HS = "ieu-b-14")
N_OUT    <- c(EPI_ALL = 44889, EPI_FOCAL = 39926, EPI_HS = 30480)

## ---- 1. 读取 MetaBrain 文件（自动适配有无 MetaBeta/MetaSE）----
read_metabrain <- function(f) {
  d <- fread(f, sep = "\t")
  d <- d[toupper(d$GeneSymbol) == toupper(GENE_SYMBOL), ]
  if (!nrow(d)) stop("文件中未找到 ", GENE_SYMBOL)
  ## SNP 格式: 14:54673482:rs2057370:C_T
  f4  <- tstrsplit(d$SNP, ":")
  a12 <- tstrsplit(f4[[4]], "_")
  has_beta <- all(c("MetaBeta", "MetaSE") %in% names(d))
  data.frame(
    SNP = f4[[3]],
    effect_allele = a12[[1]],
    other_allele  = a12[[2]],
    eaf  = d$SNPEffectAlleleFreq,
    beta = if (has_beta) d$MetaBeta else d$MetaPZ,
    se   = if (has_beta) d$MetaSE   else 1,
    pval = d$MetaP,
    gene = GENE_SYMBOL, id.exposure = NA, exposure = NA,
    stringsAsFactors = FALSE
  )
}

load_exposure <- function(name, path) {
  if (!file.exists(path)) { cat("[跳过] 文件不存在:", path, "\n"); return(NULL) }
  cat("---- 读取", name, "----\n")
  exp <- read_metabrain(path)
  cat(name, "LGALS3 cis SNP 总数:", nrow(exp),
      " 最小P:", signif(min(exp$pval), 3), "\n")
  exp$id.exposure <- name
  exp$exposure    <- paste0("LGALS3 brain eQTL (", name, ")")
  ## 筛选显著工具 + LD clumping（OpenGWAS API）
  exp <- exp[exp$pval < 1e-5, ]
  if (nrow(exp) < 2) {
    cat("[提示]", name, "显著工具 SNP 不足（<2），该来源不进入 MR\n")
    return(NULL)
  }
  exp <- tryCatch({
    cl <- ld_clump(dplyr::tibble(rsid = exp$SNP, pval = exp$pval, id = name),
                   clump_r2 = 0.1, clump_kb = 250)
    exp[exp$SNP %in% cl$rsid, ]
  }, error = function(e) {
    cat("[警告] clump 失败，使用全部 SNP:", conditionMessage(e), "\n")
    exp
  })
  cat(name, "独立工具 SNP 数:", nrow(exp), "\n")
  format_data(exp, type = "exposure")
}

exposures <- list()
for (nm in names(FILES)) {
  exposures[[nm]] <- tryCatch(
    load_exposure(nm, FILES[[nm]]),
    error = function(e) {
      cat("[失败]", nm, ":", conditionMessage(e), "\n")
      NULL
    }
  )
}
exposures <- exposures[!vapply(exposures, is.null, logical(1))]
if (!length(exposures)) stop("没有可用的脑组织暴露（可能均无显著 cis 信号）")
cat("\n可用暴露:", paste(names(exposures), collapse = ", "), "\n\n")

## ---- 2. 逐暴露 × 逐结局跑 MR ----
run_mr <- function(exp_dat, out_id, out_label, n_out) {
  out <- tryCatch(extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id),
                  error = function(e) NULL)
  if (is.null(out) || !nrow(out)) {
    cat("  [跳过] 结局无匹配 SNP:", out_id, "\n"); return(NULL)
  }
  dat <- harmonise_data(exp_dat, out)
  dat <- dat[dat$mr_keep, ]
  if (nrow(dat) < 1) return(NULL)
  res <- if (nrow(dat) >= 3) {
    mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  } else {
    mr(dat, method_list = "mr_wald_ratio")
  }
  res$outcome_n     <- n_out
  res$outcome_label <- out_label
  res$n_snp         <- nrow(dat)
  res
}

all_res <- list()
for (nm in names(exposures)) {
  for (oc in names(OUTCOMES)) {
    cat("====", nm, "->", oc, "====\n")
    key <- paste(nm, oc, sep = "_")
    all_res[[key]] <- tryCatch(
      run_mr(exposures[[nm]], OUTCOMES[oc], oc, N_OUT[oc]),
      error = function(e) { cat("  [失败]", conditionMessage(e), "\n"); NULL })
  }
}

res_df <- bind_rows(all_res)
write.csv(res_df, file.path(OUT_DIR, "16_brain_eQTL_MR_results.csv"), row.names = FALSE)

## ---- 3. 关键结果汇总（IVW / Wald）----
cat("\n========== 脑组织 eQTL MR 关键结果 ==========\n")
key <- res_df[res_df$method %in% c("Inverse variance weighted", "Wald ratio"), ]
print(key[, c("exposure", "outcome_label", "method", "nsnp", "b", "se", "pval")])

## ---- 4. 与血浆 pQTL 结果并排对照 ----
plasma <- tryCatch(read.csv(file.path(OUT_DIR, "15_MR_LGALS3_EPI_HS_vs_ieu-b-14.csv")),
                   error = function(e) NULL)
if (!is.null(plasma)) {
  pv <- plasma[plasma$method == "Inverse variance weighted", ]
  cat("\n[对照] 血浆 pQTL → HS（脚本15）: b =", round(pv$b[1], 5),
      " P =", signif(pv$pval[1], 3), "\n")
}

cat("\n结果已保存:", file.path(OUT_DIR, "16_brain_eQTL_MR_results.csv"), "\n")
cat("\n[解读提示]\n")
cat("  皮层 → HS 阳性且血浆同向 → 血浆-脑外推成立，MR 证据升级\n")
cat("  皮层 → HS 阴性 + 海马无工具 → 脑 mRNA 层面无信号，与转录后调控结论一致\n")
cat("  皮层对 EPI_ALL/EPI_FOCAL 的模式 → 检验 HS 亚型特异性是否在脑组织复现\n")
