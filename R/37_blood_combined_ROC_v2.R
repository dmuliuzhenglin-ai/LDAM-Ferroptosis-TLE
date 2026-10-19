# ============================================================================
# Script      : 37_blood_combined_ROC_v2.R
# Purpose     : Fixed version of 37_blood_combined_ROC.R (stratified folds, explicit glm data frames, stable replicate types). Supersedes 37.
# Manuscript  : Supplementary analyses (blood ROC; Results 3.2 related)
# Input       : as 37_blood_combined_ROC.R
# Output      : results/37_blood_combined_ROC/
# Dependencies: ggplot2, pROC, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 37_blood_combined_ROC.R  【完整修复版 —— 重启 R 后整体 source】
## 目的：血液双 signature 联合建模（GSE143272, 34 未用药 TLE vs 51 健康）
##   组件 A：原四基因 signature（LPCAT2, S100A10, CD9, SYNGR2）
##   组件 B：33 号检出的血液外泌体 hub 基因（ITGAM, ITGB2, CD9, YBX1）
##   问题：联合后 AUC 是否提升？（DeLong 检验 + 分层 10 折×10 次 CV）
## 修复记录：
##   - make_folds 分层分折（避免单类折）
##   - glm 用显式 data.frame(y=, s=)（避免 newdata 行数不符）
##   - n_rep 不遮蔽 rep；replicate + vapply 固定返回类型
##   - score_list 4 项与 rocs 按位置一一对应（含 NULL 占位）
## 输出：results/37_blood_combined_ROC/
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse); library(pROC); library(ggplot2)
})
set.seed(2024)

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/37_blood_combined_ROC", showWarnings = FALSE, recursive = TRUE)

## ---------- 1. 数据（与 33 号完全一致） ----------
f_e <- list.files("results", pattern = "08_GSE143272_normalized\\.csv$",
                  recursive = TRUE, full.names = TRUE)[1]
stopifnot(!is.na(f_e))
e <- read_csv(f_e, show_col_types = FALSE)
mat <- as.matrix(e[, -1]); rownames(mat) <- e[[1]]
if (anyDuplicated(rownames(mat))) {
  mat <- as.data.frame(mat) %>% rownames_to_column("g") %>%
    group_by(g) %>% summarise(across(everything(), mean), .groups = "drop") %>%
    column_to_rownames("g") %>% as.matrix()
}
f_m <- list.files("results", pattern = "08_GSE143272.*(meta|sample|group|pheno)",
                  recursive = TRUE, full.names = TRUE)[1]
stopifnot(!is.na(f_m))
meta <- read_csv(f_m, show_col_types = FALSE)
scol <- names(meta)[sapply(meta, function(x) sum(as.character(x) %in% colnames(mat)) > 30)][1]
gcol <- names(meta)[sapply(meta, function(x)
  any(grepl("naive|treated|healthy|control", as.character(x), ignore.case = TRUE)))][1]
grp <- set_names(as.character(meta[[gcol]]), as.character(meta[[scol]]))
## 注意：原数据标签拼写为 "Drug-nave"（漏 i），grepl("naive") 抓不到 → 用反向逻辑
grp[] <- dplyr::case_when(
  grepl("healthy|control", grp, ignore.case = TRUE) ~ "Healthy",
  grepl("treated|carbamazepine|cbz", grp, ignore.case = TRUE) ~ NA_character_,
  TRUE ~ "Naive")
keep <- names(grp)[!is.na(grp) & names(grp) %in% colnames(mat)]
grp  <- factor(grp[keep], levels = c("Healthy","Naive"))   # Naive=阳性
mat  <- mat[, keep, drop = FALSE]
message("分析样本: ", sum(grp=="Naive"), " Naive vs ", sum(grp=="Healthy"), " Healthy")
stopifnot(sum(grp=="Naive") == 34, sum(grp=="Healthy") == 51)

## ---------- 2. 特征 ----------
sig_orig  <- c("LPCAT2","S100A10","CD9","SYNGR2")          # 组件A（原signature）
sig_blood <- c("ITGAM","ITGB2","CD9","YBX1")               # 组件B（33号hub基因）
sig_all   <- union(sig_orig, sig_blood)                    # 7 基因（CD9去重）
snare_genes <- c("SNAP23","STX4","STX6","STX7","VAMP3","VAMP7","VTI1B")
message("各特征检出: 原sig=", sum(sig_orig %in% rownames(mat)),
        "/4, 血液=", sum(sig_blood %in% rownames(mat)),
        "/4, 联合=", sum(sig_all %in% rownames(mat)), "/7")

miss <- setdiff(sig_all, rownames(mat))
if (length(miss)) warning("缺失基因: ", paste(miss, collapse = ", "))

## z-score 标准化（按基因跨样本）
zmat <- t(scale(t(mat[rownames(mat) %in% c(sig_all, snare_genes), , drop = FALSE])))
zmat[is.na(zmat)] <- 0

score_orig  <- colMeans(zmat[intersect(sig_orig, rownames(zmat)), , drop = FALSE])
score_blood <- colMeans(zmat[intersect(sig_blood, rownames(zmat)), , drop = FALSE])
score_snare <- colMeans(zmat[intersect(snare_genes, rownames(zmat)), , drop = FALSE])

## ---------- 3. 模型 ----------
y <- as.integer(grp == "Naive")
X <- as.data.frame(t(zmat[intersect(sig_all, rownames(zmat)), , drop = FALSE]))  # 样本 × 7

fit_all <- glm(y ~ ., data = cbind(y = y, X), family = binomial)

## 分层 k 折（每折保证两类都有）
make_folds <- function(y, k = 10){
  f <- integer(length(y))
  for (cl in unique(y)) {
    idx <- which(y == cl)
    f[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  f
}

## 单分数模型 CV（10 折 × 10 次重复）
cv_auc <- function(score, y, k = 10, n_rep = 10){
  aucs <- replicate(n_rep, {
    f <- make_folds(y, k)
    mean(vapply(seq_len(k), function(i){
      tr <- f != i; te <- f == i
      p <- predict(glm(y ~ s, data = data.frame(y = y[tr], s = score[tr]), family = binomial),
                   newdata = data.frame(s = score[te]), type = "response")
      as.numeric(pROC::auc(pROC::roc(y[te], p, quiet = TRUE)))
    }, numeric(1)))
  })
  c(CV_AUC = mean(aucs), CV_sd = sd(aucs))
}

## 多基因 logit 模型 CV
cv_auc_logit <- function(X, y, k = 10, n_rep = 10){
  aucs <- replicate(n_rep, {
    f <- make_folds(y, k)
    mean(vapply(seq_len(k), function(i){
      tr <- f != i; te <- f == i
      p <- predict(glm(y ~ ., data = cbind(y = y[tr], X[tr, , drop = FALSE]), family = binomial),
                   newdata = X[te, , drop = FALSE], type = "response")
      as.numeric(pROC::auc(pROC::roc(y[te], p, quiet = TRUE)))
    }, numeric(1)))
  })
  c(CV_AUC = mean(aucs), CV_sd = sd(aucs))
}

rocs <- list(
  `原signature_4基因`   = roc(y, score_orig,  quiet = TRUE, direction = "<"),
  `血液hub_4基因`       = roc(y, score_blood, quiet = TRUE, direction = "<"),
  `联合_7基因logit`     = roc(y, predict(fit_all, type = "response"), quiet = TRUE, direction = "<"),
  `SNARE模块评分`       = roc(y, score_snare, quiet = TRUE, direction = "<")
)
## ---------- 汇总表（for 循环版，已验证可用） ----------
n_mod <- length(rocs)
v_model <- names(rocs)
v_auc <- v_ci_lo <- v_ci_hi <- v_cv <- v_cvsd <- numeric(n_mod)

for (i in seq_len(n_mod)) {
  a  <- rocs[[i]]
  ci <- as.numeric(pROC::ci.auc(a))
  cv <- if (i == 3) cv_auc_logit(X, y) else cv_auc(score_list[[i]], y)
  v_auc[i]  <- as.numeric(pROC::auc(a))
  v_ci_lo[i] <- ci[1]; v_ci_hi[i] <- ci[3]
  v_cv[i]   <- unname(cv["CV_AUC"]); v_cvsd[i] <- unname(cv["CV_sd"])
}

tab <- data.frame(model = v_model, AUC = v_auc, CI_lo = v_ci_lo,
                  CI_hi = v_ci_hi, CV_AUC = v_cv, CV_sd = v_cvsd)
print(tab)
write.csv(tab, "results/37_blood_combined_ROC/37_roc_models.csv", row.names = FALSE)

## DeLong 三组对比
dl31 <- pROC::roc.test(rocs[[3]], rocs[[1]], method = "delong")
dl21 <- pROC::roc.test(rocs[[2]], rocs[[1]], method = "delong")
dl32 <- pROC::roc.test(rocs[[3]], rocs[[2]], method = "delong")
cat("DeLong 联合 vs 原signature: P =", dl31$p.value, "\n")
cat("DeLong 血液hub vs 原signature: P =", dl21$p.value, "\n")
cat("DeLong 联合 vs 血液hub: P =", dl32$p.value, "\n")
## 与 rocs 按位置一一对应；i=3 走 cv_auc_logit，占位 NULL
score_list <- list(score_orig, score_blood, NULL, score_snare)

tab <- bind_rows(lapply(seq_along(names(rocs)), function(i){
  nm <- names(rocs)[i]
  a  <- rocs[[i]]
  ci <- as.numeric(pROC::ci.auc(a))
  cv <- if (i == 3) cv_auc_logit(X, y) else cv_auc(score_list[[i]], y)
  data.frame(model = nm,
             AUC = as.numeric(a),
             CI_lo = ci[1], CI_hi = ci[3],
             CV_AUC = unname(cv["CV_AUC"]), CV_sd = unname(cv["CV_sd"]),
             row.names = NULL)
}))
write_csv(tab, "results/37_blood_combined_ROC/37_roc_models.csv")
print(tab)

## DeLong：联合 vs 原 signature；hub vs 原 signature
dl31 <- roc.test(rocs[[3]], rocs[[1]], method = "delong")
dl21 <- roc.test(rocs[[2]], rocs[[1]], method = "delong")
dl32 <- roc.test(rocs[[3]], rocs[[2]], method = "delong")
message(sprintf("DeLong 联合 vs 原signature: P = %.3g", dl31$p.value))
message(sprintf("DeLong 血液hub vs 原signature: P = %.3g", dl21$p.value))
message(sprintf("DeLong 联合 vs 血液hub: P = %.3g", dl32$p.value))

## ---------- 4. ROC 图 ----------
cols <- c("#8C98A8", "#1783FE", "#D64545", "#67B8A6")
pdf("results/37_blood_combined_ROC/37_roc_curves.pdf", width = 6.5, height = 6)
pROC::plot.roc(rocs[[1]], col = cols[1], lwd = 2, legacy.axes = TRUE,
               main = "GSE143272 blood: Naive-TLE vs Healthy")
for (i in 2:4) pROC::plot.roc(rocs[[i]], col = cols[i], lwd = ifelse(i==3, 2.5, 2), add = TRUE)
legend("bottomright",
       legend = paste0(names(rocs), " (AUC=",
                       round(vapply(rocs, function(r) as.numeric(pROC::auc(r)), numeric(1)), 3), ")"),
       col = cols, lwd = 2, cex = 0.8, bty = "n")
dev.off()

## ---------- 5. 自动结论 ----------
best <- tab %>% slice_max(CV_AUC, n = 1, with_ties = FALSE)
verdict <- if (dl31$p.value < 0.05 && tab$CV_AUC[3] > tab$CV_AUC[1]) {
  sprintf("【联合有效】7 基因联合模型 AUC=%.2f (CV %.2f±%.2f)，显著优于原 signature (AUC=%.2f, DeLong P=%.3g)\n→ 建议实验6的血液面板采用联合 7 基因 qPCR + Gal-3 ELISA",
          tab$AUC[3], tab$CV_AUC[3], tab$CV_sd[3], tab$AUC[1], dl31$p.value)
} else if (dl32$p.value >= 0.05 && tab$CV_AUC[2] >= tab$CV_AUC[3]) {
  sprintf("【hub 面板足够】血液hub 4 基因 CV AUC=%.2f 不低于联合模型（%.2f），且 DeLong 联合 vs hub P=%.3g 不显著\n→ 成本最优：实验6面板 = hub 4 基因 + LGALS3 共 5 个指标",
          tab$CV_AUC[2], tab$CV_AUC[3], dl32$p.value)
} else {
  sprintf("【联合表观提升但需验证】联合 AUC=%.2f vs 原 %.2f, DeLong P=%.3g; CV AUC=%.2f\n→ 写 nominal 增益，实验6两个面板都测，用真实数据再选",
          tab$AUC[3], tab$AUC[1], dl31$p.value, tab$CV_AUC[3])
}
txt <- c("==== 分析 37 结论 ====", verdict, "",
         "表观 AUC（乐观）仅供参考，决策以 CV AUC 为准", "",
         "各模型:",
         paste(sprintf("%s: AUC=%.3f (%.3f-%.3f), CV=%.3f±%.3f",
                       tab$model, tab$AUC, tab$CI_lo, tab$CI_hi, tab$CV_AUC, tab$CV_sd),
               collapse = "\n"), "",
         sprintf("DeLong 联合 vs 原signature: P = %.3g", dl31$p.value),
         sprintf("DeLong 血液hub vs 原signature: P = %.3g", dl21$p.value),
         sprintf("DeLong 联合 vs 血液hub: P = %.3g", dl32$p.value))
cat(txt, file = "results/37_blood_combined_ROC/37_summary.txt", sep = "\n")
cat(txt, sep = "\n")
message("完成 → results/37_blood_combined_ROC/")
## ---------- DeLong 三组对比 ----------
dl31 <- pROC::roc.test(rocs[[3]], rocs[[1]], method = "delong")
dl21 <- pROC::roc.test(rocs[[2]], rocs[[1]], method = "delong")
dl32 <- pROC::roc.test(rocs[[3]], rocs[[2]], method = "delong")
cat("DeLong 联合 vs 原signature: P =", dl31$p.value, "\n")
cat("DeLong 血液hub vs 原signature: P =", dl21$p.value, "\n")
cat("DeLong 联合 vs 血液hub: P =", dl32$p.value, "\n")

## ---------- ROC 曲线图 ----------
cols <- c("#8C98A8", "#1783FE", "#D64545", "#67B8A6")
pdf("results/37_blood_combined_ROC/37_roc_curves.pdf", width = 6.5, height = 6)
pROC::plot.roc(rocs[[1]], col = cols[1], lwd = 2, legacy.axes = TRUE,
               main = "GSE143272 blood: Naive-TLE vs Healthy")
for (i in 2:4) pROC::plot.roc(rocs[[i]], col = cols[i], lwd = ifelse(i==3, 2.5, 2), add = TRUE)
legend("bottomright",
       legend = paste0(names(rocs), " (AUC=",
                       round(vapply(rocs, function(r) as.numeric(pROC::auc(r)), numeric(1)), 3), ")"),
       col = cols, lwd = 2, cex = 0.8, bty = "n")
dev.off()

## ---------- 写 summary ----------
txt <- c("==== 分析 37 结论 ====", "",
         paste0("各模型（表观 AUC + 10折x10次 CV）:"),
         paste(sprintf("%s: AUC=%.3f (%.3f-%.3f), CV=%.3f±%.3f",
                       tab$model, tab$AUC, tab$CI_lo, tab$CI_hi, tab$CV_AUC, tab$CV_sd),
               collapse = "\n"), "",
         sprintf("DeLong 联合 vs 原signature: P = %.4g", dl31$p.value),
         sprintf("DeLong 血液hub vs 原signature: P = %.4g", dl21$p.value),
         sprintf("DeLong 联合 vs 血液hub: P = %.4g", dl32$p.value))
cat(txt, file = "results/37_blood_combined_ROC/37_summary.txt", sep = "\n")
cat(txt, sep = "\n")
