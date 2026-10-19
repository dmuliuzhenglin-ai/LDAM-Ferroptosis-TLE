# ============================================================================
# Script      : 04_machine_learning_hub_ROC.R
# Purpose     : Hub-gene screening by LASSO + SVM-RFE + random forest intersection; ROC in discovery and external blood cohorts.
# Manuscript  : Methods 2.3 / Results 3.2
# Input       : DE results of 01/02, WGCNA modules of 03, genesets, GSE143272 matrix (from 08 Part A)
# Output      : results/04_* hub genes, model objects, ROC tables/figures
# Dependencies: caret, data.table, dplyr, e1071, ggplot2, glmnet, pROC, randomForest
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

# ==============================================================================
# 04_machine_learning_hub_ROC.R
# ① DE-LDGs：两个 bulk 队列 DEG ∩ LDAM 基因集 ∩ WGCNA 关键模块
# ② LASSO + SVM-RFE + 随机森林三算法交集筛 hub 基因
# ③ ROC：发现队列 GSE168375 + 外部验证（GSE143272 外周血：未用药癫痫 vs 健康）
# 依赖：先运行 00 / 01 / 02 / 03，以及 08 中的 GSE143272 预处理部分
# ==============================================================================

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results"); FIG_DIR <- file.path(BASE_DIR, "figures")
GS_DIR   <- file.path(BASE_DIR, "genesets")

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(glmnet); library(e1071)
  library(randomForest); library(pROC); library(ggplot2); library(caret)
})
set.seed(2026)

LDAM_genes <- readLines(file.path(GS_DIR, "LDAM_geneset.txt"))

## ---------- 4.1 取交集得到候选基因（DE-LDGs ∩ WGCNA 模块，方向一致性过滤） ----------
de1 <- fread(file.path(OUT_DIR, "01_GSE134697_DEG_Hippocampus_vs_Neocortex.csv"))
deg1 <- de1 |> filter(padj < 0.05 & abs(log2FoldChange) > 1) |> pull(gene)

# 主交集：GSE134697 配对 DEGs ∩ LDAM 基因集（配对设计功效最高，作为主筛选）
DE_LDGs <- intersect(deg1, LDAM_genes)

# 验证队列方向一致性过滤：在 GSE190451（组织级）或 GSE168375（星胶）中
# 方向一致且名义显著的基因优先保留；两个队列都缺失/不一致的基因剔除
val_files <- c(GSE190451 = file.path(OUT_DIR, "02_GSE190451_DEG_TLE_vs_Control.csv"),
               GSE168375 = file.path(OUT_DIR, "02_GSE168375_DEG_Epilepsy_vs_Control.csv"))
consistent <- DE_LDGs
for (vf in val_files[file.exists(val_files)]) {
  dv <- fread(vf)
  dv_sub <- dv |> filter(gene %in% DE_LDGs) |>
    select(gene, log2FoldChange, pval = pvalue)
  de1_sub <- de1 |> filter(gene %in% DE_LDGs) |>
    select(gene, main_lfc = log2FoldChange)
  chk <- de1_sub |> inner_join(dv_sub, by = "gene") |>
    mutate(consistent = sign(main_lfc) == sign(log2FoldChange) & pval < 0.05)
  pass <- chk |> filter(consistent) |> pull(gene)
  message(basename(vf), " 方向一致基因数: ", length(pass), "/", nrow(chk))
  if (length(pass) >= 5) consistent <- intersect(consistent, pass)
}
DE_LDGs <- consistent

wgcna_files <- list.files(OUT_DIR, pattern = "03_WGCNA_top_module_.*_genes.txt",
                          full.names = TRUE)
module_genes <- if (length(wgcna_files) > 0) readLines(wgcna_files[1]) else character(0)

if (length(module_genes) > 0) {
  candidates <- intersect(DE_LDGs, module_genes)
  if (length(candidates) < 3) {           # 交集过小时回退，避免流程中断
    message("与 WGCNA 模块交集过少，回退为 DE-LDGs")
    candidates <- DE_LDGs
  }
} else candidates <- DE_LDGs
message("候选基因数: ", length(candidates))
writeLines(candidates, file.path(OUT_DIR, "04_candidate_genes.txt"))

## ---------- 4.2 机器学习训练集：GSE168375（31 癫痫 vs 12 对照） ----------
expr_vst <- fread(file.path(OUT_DIR, "02_GSE168375_vst_normalized.csv")) |>
  tibble::column_to_rownames("gene") |> as.matrix()
meta <- fread(file.path(OUT_DIR, "02_GSE168375_sample_meta.csv"))

use_genes <- intersect(candidates, rownames(expr_vst))
X <- t(expr_vst[use_genes, meta$sample])
y <- factor(meta$group, levels = c("Control", "Epilepsy"))

## ---------- 4.3 LASSO ----------
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5)
pdf(file.path(FIG_DIR, "04_LASSO_cv.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2)); plot(cvfit); plot(cvfit$glmnet.fit, xvar = "lambda"); dev.off()
lasso_genes <- rownames(coef(cvfit, s = "lambda.min"))[
  which(coef(cvfit, s = "lambda.min")[, 1] != 0)][-1]
message("LASSO 选中: ", length(lasso_genes), " 个基因")

## ---------- 4.4 SVM-RFE ----------
svm_rfe <- rfe(x = X, y = y, sizes = c(1:min(20, ncol(X))),
               rfeControl = rfeControl(functions = caretFuncs, method = "cv", number = 5),
               methods = "svmRadial")
svm_genes <- predictors(svm_rfe)
message("SVM-RFE 选中: ", length(svm_genes), " 个基因")

## ---------- 4.5 随机森林 ----------
rf <- randomForest(x = X, y = y, ntree = 1000, importance = TRUE)
imp <- importance(rf, type = 2)     # MeanDecreaseGini
rf_genes <- rownames(imp)[imp[, 1] > quantile(imp[, 1], 0.5)]   # 重要性前 50%
message("随机森林选中: ", length(rf_genes), " 个基因")

## ---------- 4.6 三算法交集 = hub 基因 ----------
hub_genes <- Reduce(intersect, list(lasso_genes, svm_genes, rf_genes))
if (length(hub_genes) == 0) {  # 严格交集为空时取"至少命中 2/3"的基因
  tab <- table(c(lasso_genes, svm_genes, rf_genes))
  hub_genes <- names(tab)[tab >= 2]
  message("三算法交集为空，回退为 ≥2/3 算法命中")
}
message("Hub 基因: ", paste(hub_genes, collapse = ", "))
writeLines(hub_genes, file.path(OUT_DIR, "04_hub_genes.txt"))
ml_lists <- list(LASSO = lasso_genes, SVM_RFE = svm_genes, RF = rf_genes)
fwrite(rbindlist(lapply(names(ml_lists),
                        function(m) data.frame(method = m, gene = ml_lists[[m]])),
                 fill = TRUE),
       file.path(OUT_DIR, "04_ML_selected_genes.csv"))

## ---------- 4.7 ROC 评估 ----------
# 辅助函数：hub 基因逻辑回归模型 → ROC
roc_eval <- function(expr_mat, sample_meta, genes, title, outfile) {
  g <- intersect(genes, rownames(expr_mat))
  Xv <- t(expr_mat[g, sample_meta$sample, drop = FALSE])
  df <- data.frame(group = factor(sample_meta$group, levels = c("Control", "Epilepsy")), Xv)
  fit <- glm(group ~ ., data = df, family = binomial)
  prob <- predict(fit, type = "response")
  roc_obj <- roc(df$group, prob, quiet = TRUE)
  p <- ggroc(roc_obj, color = "#2C5F8A", linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 1, linetype = 2, color = "grey60") +
    annotate("text", x = 0.6, y = 0.2,
             label = paste0("AUC = ", round(auc(roc_obj), 3),
                            " (95%CI ", round(ci.auc(roc_obj)[1], 3), "-",
                            round(ci.auc(roc_obj)[3], 3), ")")) +
    labs(title = title) + theme_bw()
  ggsave(outfile, p, width = 6, height = 6)
  return(auc(roc_obj))
}

# 发现队列 ROC
auc_train <- roc_eval(expr_vst, meta, hub_genes,
                      "GSE168375 训练队列 ROC（星形胶质细胞）",
                      file.path(FIG_DIR, "04_ROC_discovery_GSE168375.pdf"))
message("训练队列 AUC = ", round(auc_train, 3))

# 外部验证 1：GSE190451 组织级验证队列（TLE vs 对照，3v3）
f451 <- file.path(OUT_DIR, "02_GSE190451_normalized.csv")
if (file.exists(f451)) {
  expr_451 <- fread(f451) |> tibble::column_to_rownames("gene") |> as.matrix()
  meta_451 <- fread(file.path(OUT_DIR, "02_GSE190451_sample_meta.csv"))
  meta_451$group <- ifelse(meta_451$group == "TLE", "Epilepsy", "Control")
  auc_451 <- roc_eval(expr_451, meta_451, hub_genes,
                      "GSE190451 组织级验证队列 ROC",
                      file.path(FIG_DIR, "04_ROC_validation_GSE190451.pdf"))
  message("GSE190451 验证 AUC = ", round(auc_451, 3))
} else message("未找到 GSE190451 矩阵，请先运行 02 脚本 Part A")

# 外部验证 2：GSE143272 外周血（未用药癫痫 vs 健康对照；需先运行 08 脚本 Part A）
peri_file <- file.path(OUT_DIR, "08_GSE143272_normalized.csv")
peri_meta_file <- file.path(OUT_DIR, "08_GSE143272_sample_meta.csv")
if (file.exists(peri_file)) {
  expr_peri <- fread(peri_file) |> tibble::column_to_rownames("gene") |> as.matrix()
  meta_peri <- fread(peri_meta_file)
  auc_ext <- roc_eval(expr_peri, meta_peri, hub_genes,
                      "GSE143272 外周血验证队列 ROC（未用药癫痫 vs 健康）",
                      file.path(FIG_DIR, "04_ROC_validation_GSE143272_blood.pdf"))
  message("外周血验证队列 AUC = ", round(auc_ext, 3))
} else {
  message("未找到 GSE143272 预处理文件，请先运行 08 脚本后再执行本段 ROC 验证")
}

message("=== 04_machine_learning_hub_ROC.R 运行完毕 ===")
