# ============================================================================
# Script      : 37_blood_combined_ROC.R
# Purpose     : Blood dual-signature joint modelling (original 4-gene signature + blood exosome hub genes): apparent and repeated 10-fold CV AUC, DeLong test.
# Manuscript  : Supplementary analyses (blood ROC; Results 3.2 related)
# Input       : GSE143272 matrix; gene lists from 04/33
# Output      : results/37_blood_combined_ROC/
# Dependencies: ggplot2, pROC, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 37_blood_combined_ROC.R
## 目的：血液双 signature 联合建模（GSE143272, 34 未用药 TLE vs 51 健康）
##   组件 A：原四基因 signature（LPCAT2, S100A10, CD9, SYNGR2）—— 主文 AUC=0.74
##   组件 B：33 号检出的血液外泌体 hub 基因（ITGAM, ITGB2, CD9, YBX1）
##   问题：联合后 AUC 是否提升？（DeLong 检验 + 交叉验证）
## 设计要点：
##   1. 报告【表观 AUC】（optimistic）和【10 折×10 次重复 CV AUC】（honest）两套
##   2. 主要对比：原 signature vs 联合 7 基因（CD9 去重）
##   3. 次要对比：组件 B 单独、SNARE_fusion 模块评分单独
## 输出：results/37_blood_combined_ROC/
##   37_roc_models.csv      各模型 AUC + CI + CV AUC
##   37_roc_curves.pdf      ROC 曲线
##   37_summary.txt         自动结论（含实验6预实验建议）
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
grp[] <- dplyr::case_when(
  grepl("naive", grp, ignore.case = TRUE) ~ "Naive",
  grepl("healthy|control", grp, ignore.case = TRUE) ~ "Healthy",
  TRUE ~ NA_character_)
keep <- names(grp)[!is.na(grp) & grp %in% c("Naive","Healthy") & names(grp) %in% colnames(mat)]
grp  <- factor(grp[keep], levels = c("Healthy","Naive"))   # Naive=阳性
mat  <- mat[, keep, drop = FALSE]
message("分析样本: ", sum(grp=="Naive"), " Naive vs ", sum(grp=="Healthy"), " Healthy")

## ---------- 2. 特征 ----------
sig_orig  <- c("LPCAT2","S100A10","CD9","SYNGR2")          # 组件A（原signature）
sig_blood <- c("ITGAM","ITGB2","CD9","YBX1")               # 组件B（33号hub基因）
sig_all   <- union(sig_orig, sig_blood)                    # 7 基因（CD9去重）
snare_genes <- c("SNAP23","STX4","STX6","STX7","VAMP3","VAMP7","VTI1B")
message("各特征检出: 原sig=", sum(sig_orig %in% rownames(mat)),
        "/4, 血液=", sum(sig_blood %in% rownames(mat)),
        "/4, 联合=", length(sig_all), "基因中", sum(sig_all %in% rownames(mat)), "检出")

## 缺失基因告警但不中断
miss <- setdiff(sig_all, rownames(mat))
if (length(miss)) warning("缺失基因: ", paste(miss, collapse = ", "))

## z-score 标准化（按基因跨样本）
zmat <- t(scale(t(mat[rownames(mat) %in% sig_all | rownames(mat) %in% snare_genes, , drop = FALSE])))
zmat[is.na(zmat)] <- 0

## 单特征分数
score_orig  <- colMeans(zmat[intersect(sig_orig, rownames(zmat)), , drop = FALSE])
score_blood <- colMeans(zmat[intersect(sig_blood, rownames(zmat)), , drop = FALSE])
score_snare <- colMeans(zmat[intersect(snare_genes, rownames(zmat)), , drop = FALSE])

## ---------- 3. 模型 ----------
## 逻辑回归（注意 glm 的方向：以 Naive 为阳性，用 y = (grp=="Naive")）
y <- as.integer(grp == "Naive")
X <- t(zmat[sig_all, , drop = FALSE])   # 样本 × 7 基因

fit_orig  <- glm(y ~ score_orig,  family = binomial)
fit_blood <- glm(y ~ score_blood, family = binomial)
fit_all   <- glm(y ~ ., data = as.data.frame(X), family = binomial)

## 10 折 × 10 次重复交叉验证 AUC
cv_auc <- function(score, y, k = 10, rep = 10){
  aucs <- sapply(seq_len(rep), function(r){
    f <- sample(rep(seq_len(k), length.out = length(y)))
    mean(sapply(seq_len(k), function(i){
      tr <- f != i; te <- f == i
      p <- predict(glm(y[tr] ~ score[tr], family = binomial), newdata = data.frame(score = score[te]), type = "response")
      as.numeric(pROC::auc(pROC::roc(y[te], p, quiet = TRUE)))
    }))
  })
  c(mean = mean(aucs), sd = sd(aucs))
}
cv_auc_logit <- function(X, y, k = 10, rep = 10){
  aucs <- sapply(seq_len(rep), function(r){
    f <- sample(rep(seq_len(k), length.out = length(y)))
    mean(sapply(seq_len(k), function(i){
      tr <- f != i; te <- f == i
      p <- predict(glm(y[tr] ~ ., data = as.data.frame(X)[tr, , drop = FALSE], family = binomial),
                   newdata = as.data.frame(X)[te, , drop = FALSE], type = "response")
      as.numeric(pROC::auc(pROC::roc(y[te], p, quiet = TRUE)))
    }))
  })
  c(mean = mean(aucs), sd = sd(aucs))
}

rocs <- list(
  `原signature_4基因`   = roc(y, score_orig,  quiet = TRUE, direction = "<"),
  `血液hub_4基因`       = roc(y, score_blood, quiet = TRUE, direction = "<"),
  `联合_7基因logit`     = roc(y, predict(fit_all, type = "response"), quiet = TRUE, direction = "<"),
  `SNARE模块评分`       = roc(y, score_snare, quiet = TRUE, direction = "<")
)
tab <- bind_rows(lapply(names(rocs), function(nm){
  a <- rocs[[nm]]
  ci <- pROC::ci.auc(a)
  cv <- if (nm == "联合_7基因logit") cv_auc_logit(X, y) else cv_auc(scores <- switch(nm,
    `原signature_4基因` = score_orig, `血液hub_4基因` = score_blood, `SNARE模块评分` = score_snare), y)
  tibble(model = nm, AUC = as.numeric(a), CI_low = ci[1], CI_high = ci[3],
         CV_AUC = cv["mean"], CV_sd = cv["sd"])
}))
write_csv(tab, "results/37_blood_combined_ROC/37_roc_models.csv")
print(tab)

## DeLong：联合 vs 原 signature
dl <- roc.test(rocs[["联合_7基因logit"]], rocs[["原signature_4基因"]], method = "delong")
message(sprintf("DeLong 联合 vs 原signature: P = %.3g", dl$p.value))

## ---------- 4. ROC 图 ----------
ggroc(rocs, legacy.axes = TRUE, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  scale_color_manual(values = c("#0E3A6B","#2E86C8","#C0392B","#7FB069"),
    labels = paste0(tab$model, sprintf(" (AUC=%.2f)", tab$AUC))) +
  labs(title = "GSE143272 血液诊断模型：Naive TLE vs Healthy",
       x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  theme_classic(base_size = 12)
ggsave("results/37_blood_combined_ROC/37_roc_curves.pdf", width = 7, height = 6)

## ---------- 5. 自动结论 ----------
best <- tab %>% slice_max(AUC, n = 1, with_ties = FALSE)
dl_p <- dl$p.value
verdict <- if (dl_p < 0.05 && best$model == "联合_7基因logit" && best$CV_AUC > 0.70) {
  sprintf("【联合有效】7 基因联合模型 AUC=%.2f (CV %.2f±%.2f)，显著优于原 signature (AUC=%.2f, DeLong P=%.3g)
→ 建议实验6的血液面板采用联合 7 基因 qPCR + Gal-3 ELISA", best$AUC, best$CV_AUC, best$CV_sd,
          tab$AUC[tab$model=="原signature_4基因"], dl_p)
} else if (best$model != "联合_7基因logit") {
  sprintf("【原signature仍最优】最佳模型=%s (AUC=%.2f, CV %.2f)；联合无增益 → 实验6维持原 4 基因面板",
          best$model, best$AUC, best$CV_AUC)
} else {
  sprintf("【联合表观提升但需验证】联合 AUC=%.2f vs 原 %.2f，DeLong P=%.3g；CV AUC=%.2f
→ 写 nominal 增益，实验6两个面板都测，用真实数据再选", best$AUC,
          tab$AUC[tab$model=="原signature_4基因"], dl_p, best$CV_AUC)
}
cat(c("==== 分析 37 结论 ====", verdict, "",
      "表观 AUC（乐观）仅供参考，决策以 CV AUC 为准", "",
      "各模型:", paste(sprintf("%s: AUC=%.2f (%.2f-%.2f), CV=%.2f±%.2f",
        tab$model, tab$AUC, tab$CI_low, tab$CI_high, tab$CV_AUC, tab$CV_sd), collapse = "\n")),
    file = "results/37_blood_combined_ROC/37_summary.txt", sep = "\n")
cat(readLines("results/37_blood_combined_ROC/37_summary.txt"), sep = "\n")
