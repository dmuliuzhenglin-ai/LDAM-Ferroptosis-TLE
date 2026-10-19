# ============================================================================
# Script      : 15_MR分析.R
# Purpose     : Two-sample Mendelian randomization: Galectin-3 pQTL (deCODE), LDL/TG and hub-gene eQTL vs epilepsy GWAS (OpenGWAS).
# Manuscript  : Methods 2.10 / Results 3.11
# Input       : OpenGWAS API (ieu-b-8/-10/-14 outcomes); deCODE pQTL; eQTLGen
# Output      : results/15_* MR tables
# Dependencies: TwoSampleMR, coloc, data.table, dplyr, ieugwasr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

# [SECURITY] OpenGWAS JWT redacted for publication. Provide your own token via
#   .Renviron (OPENGWAS_JWT=...) or Sys.setenv() before running. See ieugwasr docs.
if (!nzchar(Sys.getenv("OPENGWAS_JWT"))) stop("Set OPENGWAS_JWT first (see ieugwasr::get_opengwas_jwt / .Renviron)")
ieugwasr::user()
## ============================================================
## 15_MR分析.R（直连版，免全目录查询）
## 结局：ieu-b-8 全癫痫 / ieu-b-10 focal / ieu-b-14 focal+海马硬化（ILAE 2018）
## 暴露：Galectin-3 pQTL(deCODE) / LDL / TG / hub基因 eQTL(eQTLGen)
## 前置：OPENGWAS_JWT 已配置并重启 R
## ============================================================
suppressPackageStartupMessages({
  library(TwoSampleMR); library(ieugwasr); library(dplyr); library(data.table)
})
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results"); dir.create(OUT_DIR, showWarnings = FALSE)

## ---------- 0. 数据集 ID（已核实）+ 验证 ----------
IDS <- list(
  EPI_ALL   = "ieu-b-8",            # epilepsy, all documented cases (ILAE 2018)
  EPI_FOCAL = "ieu-b-10",           # focal epilepsy
  EPI_HS    = "ieu-b-14",           # focal epilepsy + hippocampal sclerosis ← HTLE 专属
  GAL3_PQTL = "prot-c-3066_12_1",   # Galectin-3, deCODE pQTL
  LDL       = "ieu-a-300",          # LDL cholesterol (GLGC)
  TG        = "ieu-a-302",          # Triglycerides (GLGC)
  SYNGR2_eQ = "eqtl-a-ENSG00000108639",  # eQTLGen blood cis-eQTL
  HSPB1_eQ  = "eqtl-a-ENSG00000106211"
)
cat("===== 0. 验证数据集 =====\n")
info <- lapply(names(IDS), function(nm){
  ii <- tryCatch(gwasinfo(IDS[[nm]]), error = function(e) NULL)
  if(is.null(ii) || nrow(ii) == 0){ cat(sprintf("  %-10s %s 【查询失败/不存在】\n", nm, IDS[[nm]])); return(NULL) }
  tr <- if(length(ii$trait) > 0) ii$trait[1] else NA
  nn <- suppressWarnings(if(length(ii$sample_size) > 0) ii$sample_size[1] else NA)
  cat(sprintf("  %-10s %s | %s | n=%s\n", nm, IDS[[nm]], tr, nn))
  data.frame(name = nm, id = IDS[[nm]], trait = ifelse(is.na(tr),"",tr),
             n = ifelse(is.na(nn), NA_integer_, as.integer(nn)))
})
info <- do.call(rbind, info[!sapply(info, is.null)])
fwrite(info, file.path(OUT_DIR, "15_MR_datasets_used.csv"))
EPI_SETS <- c(EPI_HS = IDS$EPI_HS, EPI_FOCAL = IDS$EPI_FOCAL, EPI_ALL = IDS$EPI_ALL)

run_mr <- function(exp_dat, out_id, tag){
  out <- tryCatch(extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id), error = function(e) NULL)
  if(is.null(out) || nrow(out) == 0){ cat(sprintf("  [%s vs %s] 结局无匹配 SNP\n", tag, out_id)); return(NULL) }
  dat <- harmonise_data(exp_dat, out); dat <- dat[dat$mr_keep, ]
  if(nrow(dat) == 0){ cat(sprintf("  [%s vs %s] harmonise 后为空\n", tag, out_id)); return(NULL) }
  res <- generate_odds_ratios(mr(dat, method_list = c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median")))
  cat(sprintf("--- %s vs %s: %d SNPs ---\n", tag, out_id, nrow(dat)))
  print(res[, c("method","nsnp","b","se","pval","or","or_lci95","or_uci95")])
  fwrite(res, file.path(OUT_DIR, sprintf("15_MR_%s_vs_%s.csv", tag, out_id)))
  invisible(list(dat = dat, res = res))
}

## ============================================================
## ① Galectin-3 pQTL cis-MR
## ============================================================
gal_res <- tryCatch({
  cat("\n===== ① Galectin-3 pQTL MR =====\n")
  exp_gal <- extract_instruments(outcomes = IDS$GAL3_PQTL, p1 = 5e-6, clump = TRUE, r2 = 0.1)
  cat("Galectin-3 工具变量数:", nrow(exp_gal), "\n")
  if(nrow(exp_gal) > 0){
    out_all <- lapply(names(EPI_SETS), function(epi) run_mr(exp_gal, EPI_SETS[[epi]], paste0("LGALS3_", epi)))
    names(out_all) <- names(EPI_SETS)
    out_all
  }
}, error = function(e) cat("① 出错:", conditionMessage(e), "\n"))
saveRDS(gal_res, file.path(OUT_DIR, "15_MR_LGALS3_all_outcomes.rds"))

## ============================================================
## ② 血脂 MR（LDL-C、TG）——呼应 LDAM 脂质轴
## ============================================================
tryCatch({
  cat("\n===== ② 血脂 MR =====\n")
  for(lip in c("LDL","TG")){
    exp_l <- extract_instruments(outcomes = IDS[[lip]], clump = TRUE, r2 = 0.001)
    cat(lip, "工具变量数:", nrow(exp_l), "\n")
    run_mr(exp_l, IDS$EPI_HS,    paste0(lip, "_EPI_HS"))
    run_mr(exp_l, IDS$EPI_FOCAL, paste0(lip, "_EPI_FOCAL"))
  }
}, error = function(e) cat("② 出错:", conditionMessage(e), "\n"))

## ============================================================
## ③ 共定位（LGALS3 位点区域 SNP，coloc.abf）
## ============================================================
tryCatch({
  cat("\n===== ③ 共定位 coloc（LGALS3 vs focal+HS） =====\n")
  library(coloc)
  exp_full <- extract_instruments(outcomes = IDS$GAL3_PQTL, p1 = 5e-4, clump = FALSE)
  cat("区域 SNP 数:", nrow(exp_full), "\n")
  if(nrow(exp_full) >= 20){
    out_full <- extract_outcome_data(snps = exp_full$SNP, outcomes = IDS$EPI_HS)
    dat <- harmonise_data(exp_full, out_full)
    dat <- dat[complete.cases(dat$beta.exposure, dat$beta.outcome), ]
    cat("共定位 SNP 数:", nrow(dat), "\n")
    if(nrow(dat) >= 20){
      n_exp <- info$n[info$name == "GAL3_PQTL"]; n_out <- info$n[info$name == "EPI_HS"]
      cc <- coloc.abf(
        dataset1 = list(beta = dat$beta.exposure, varbeta = dat$se.exposure^2, N = n_exp, type = "quant", snp = dat$SNP),
        dataset2 = list(beta = dat$beta.outcome,  varbeta = dat$se.outcome^2,  N = n_out, type = "cc",   snp = dat$SNP))
      print(cc$summary)
      fwrite(as.data.frame(t(cc$summary)), file.path(OUT_DIR, "15_coloc_LGALS3_HS.csv"))
      cat("解读: PP.H4 > 0.7 → 同一因果变异，MR 结果可信；PP.H3 高 → 连锁偏倚，MR 不可信\n")
    } else cat("  SNP 不足（<20），近似 coloc 做不了；需要 deCODE 全量 summary stats 离线做\n")
  } else cat("  区域 SNP 太少，coloc 跳过\n")
}, error = function(e) cat("③ 出错:", conditionMessage(e), "\n"))

## ============================================================
## ④ hub 基因 eQTL MR（eQTLGen 血液 cis-eQTL）
## ============================================================
tryCatch({
  cat("\n===== ④ hub 基因 eQTL MR =====\n")
  for(gene in c("SYNGR2","HSPB1")){
    gid <- IDS[[paste0(gene, "_eQ")]]
    exp_e <- tryCatch(extract_instruments(outcomes = gid, p1 = 5e-6, clump = TRUE, r2 = 0.1),
                      error = function(e) NULL)
    if(is.null(exp_e) || nrow(exp_e) == 0){ cat(" ", gene, "无可用 cis-eQTL 工具变量，跳过\n"); next }
    cat(gene, "工具变量数:", nrow(exp_e), "\n")
    run_mr(exp_e, IDS$EPI_HS,    paste0(gene, "_eQTL_HS"))
    run_mr(exp_e, IDS$EPI_FOCAL, paste0(gene, "_eQTL_FOCAL"))
  }
}, error = function(e) cat("④ 出错:", conditionMessage(e), "\n"))

cat("\n全部完成，结果在: ", OUT_DIR, "\n")
