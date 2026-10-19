# ============================================================================
# Script      : 07_scRNA_GSE201048_validation.R
# Purpose     : LEGACY TEMPLATE - BayesPrism deconvolution workflow carried over from an earlier glioma project; NOT part of this manuscript.
# Manuscript  : not used in manuscript (legacy file, kept for transparency)
# Input       : external glioma-project rds/TCGA files (paths parameterised, data not in this repo)
# Output      : 08_prism_theta_cox.csv (legacy)
# Dependencies: BayesPrism, Matrix, Seurat, data.table, survival
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Status      : LEGACY / DEPRECATED - kept for transparency, not part of the manuscript pipeline
# ============================================================================

## ============================================================================
## 08 BayesPrism 去卷积（最终修复版 v3）
## 输入: myeloid_final.rds + myeloid_counts_raw.rds（剔除 GSE131928 TPM 层）
##       + GSE182109_raw.rds / GSE138794_qc.rds（非髓系 counts 来源）
##       + TCGA star_counts / survival / clinical
## 输出: 08_prism_theta_cox.csv
## 用法: 重启 R 后整文件运行；检查点自动续跑
## ============================================================================

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root

options(future.globals.maxSize = 32 * 1024^3)
# [legacy] original paths pointed to an earlier glioma project; parameterised below
GEO_DIR <- Sys.getenv("GLIOMA_GEO_DIR", file.path(BASE_DIR, "data/external/glioma/GEO"))
setwd(GEO_DIR)
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
  library(BayesPrism); library(survival)
})

RDS_DIR  <- file.path(GEO_DIR, "analysis", "rds")
TCGA_DIR <- Sys.getenv("GLIOMA_TCGA_DIR", file.path(BASE_DIR, "data/external/glioma/TCGA"))
OUT_DIR  <- file.path(GEO_DIR, "analysis", "bayesprism")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ref_rds   <- file.path(OUT_DIR, "08_prism_reference.rds")
bulk_rds  <- file.path(OUT_DIR, "08_prism_bulk.rds")
prism_rds <- file.path(OUT_DIR, "08_prism_result.rds")

## 工具函数：把若干稀疏矩阵按"基因并集"对齐（缺失基因补零）后 cbind
align_union_cbind <- function(lst) {
  lst <- lst[!sapply(lst, is.null)]
  stopifnot(length(lst) > 0)
  g_all <- Reduce(union, lapply(lst, rownames))
  lst2 <- lapply(lst, function(m) {
    miss <- setdiff(g_all, rownames(m))
    if (length(miss) > 0) {
      zm <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                 dims = c(length(miss), ncol(m)),
                                 dimnames = list(miss, colnames(m)))
      m <- rbind(m, zm)
    }
    m[g_all, , drop = FALSE]
  })
  do.call(cbind, lst2)
}

## ============================================================================
## 1. 构建参考矩阵
## ============================================================================
if (file.exists(ref_rds)) {
  message("读取已保存的参考矩阵...")
  ref <- readRDS(ref_rds)
  sc_dat <- ref$sc_dat; cell_type <- ref$cell_type; cell_state <- ref$cell_state
  rm(ref); gc()
} else {

  ## ---- 1a. 标签：髓系亚型 + 非髓系大类（剔除 GSE131928 细胞）----
  myeloid <- readRDS(file.path(RDS_DIR, "myeloid_final.rds"))
  sub_lab <- setNames(as.character(myeloid$subtype_final), colnames(myeloid))
  broad_m <- setNames(as.character(myeloid$marker_based_anno), colnames(myeloid))
  rm(myeloid); gc()
  keep_mye <- !grepl("GSE131928", names(sub_lab))     # TPM 数据集，剔除
  sub_lab <- sub_lab[keep_mye]; broad_m <- broad_m[keep_mye]
  message("髓系细胞（剔除 GSE131928 后）: ", length(sub_lab))

  merged <- readRDS(file.path(RDS_DIR, "merged_annotated.rds"))
  broad_all <- setNames(as.character(merged$marker_based_anno), colnames(merged))
  rm(merged); gc()
  non_mye_lab <- broad_all[broad_all %in% c("Tumor","Astrocyte","Oligodendro",
                                            "Endothelial","Pericyte","Fibroblast",
                                            "T_cell","NK","B_cell")]
  non_mye_lab <- non_mye_lab[!grepl("GSE131928", names(non_mye_lab))]
  non_mye_lab <- non_mye_lab[setdiff(names(non_mye_lab), names(sub_lab))]
  message("非髓系细胞（剔除 GSE131928 后）: ", length(non_mye_lab))

  ## ---- 1b. 先抽样（省内存）：每 state ≤1000，Tumor ≤3000 ----
  state_all <- c(sub_lab, non_mye_lab)
  type_all  <- c(broad_m, non_mye_lab)
  set.seed(42)
  keep_cells <- unlist(lapply(split(names(state_all), state_all), function(v) {
    cap <- if (v[1] == "Tumor") 3000 else 1000
    if (length(v) > cap) sample(v, cap) else v
  }))
  message("抽样后参考细胞数: ", length(keep_cells))

  ## ---- 1c. 髓系 counts：逐层取 keep_cells，基因并集对齐 ----
  cnt_list <- readRDS(file.path(RDS_DIR, "myeloid_counts_raw.rds"))
  cnt_list <- cnt_list[!grepl("GSE131928", names(cnt_list))]   # 剔除 TPM 层
  mye_list <- lapply(cnt_list, function(m) {
    keep <- intersect(colnames(m), keep_cells)
    if (length(keep) == 0) return(NULL)
    m[, keep, drop = FALSE]
  })
  rm(cnt_list); gc()
  mye_counts <- align_union_cbind(mye_list)
  rm(mye_list); gc()
  message("髓系 counts: ", nrow(mye_counts), " 基因 x ", ncol(mye_counts), " 细胞")

  ## ---- 1d. 非髓系 counts：从原始 GSE 文件取 ----
  nm_list <- list()
  for (f in c(file.path(RDS_DIR, "GSE182109_raw.rds"),
              file.path(RDS_DIR, "GSE138794_qc.rds"))) {
    obj <- readRDS(f)
    for (ln in Layers(obj[["RNA"]], search = "^counts")) {
      m <- LayerData(obj, assay = "RNA", layer = ln)
      keep <- intersect(colnames(m), keep_cells)
      if (length(keep) == 0) next
      nm_list[[paste0(basename(f), "_", ln)]] <- m[, keep, drop = FALSE]
    }
    rm(obj); gc()
  }
  if (length(nm_list) == 0)
    stop("非髓系 counts 为空：原始文件与 merged 的细胞名可能对不上，请把上面的细胞数信息贴给我")
  nm_counts <- align_union_cbind(nm_list)
  rm(nm_list); gc()
  message("非髓系 counts: ", nrow(nm_counts), " 基因 x ", ncol(nm_counts), " 细胞")
  message("非髓系标签覆盖率: ",
          round(100 * ncol(nm_counts) / length(intersect(names(non_mye_lab), keep_cells)), 1), "%")

  ## ---- 1e. 合并髓系+非髓系（基因并集对齐）----
  sc_counts <- align_union_cbind(list(mye_counts, nm_counts))
  rm(mye_counts, nm_counts); gc()

  cells_final <- colnames(sc_counts)
  cell_state <- state_all[cells_final]
  cell_type  <- type_all[cells_final]
  names(cell_state) <- names(cell_type) <- cells_final
  stopifnot(!any(is.na(cell_state)), !any(is.na(cell_type)))
  print(table(cell_state))

  ## ---- 1f. 基因过滤 + 转置（行=细胞, 列=基因）----
  g_all <- rownames(sc_counts)
  bad <- grepl("^MT-|^RP[SL]|^MALAT1", g_all)
  det <- Matrix::rowMeans(sc_counts > 0) > 0.01
  sc_counts <- sc_counts[!bad & det, , drop = FALSE]
  message("过滤后基因数: ", nrow(sc_counts))

  sc_dat <- as.matrix(t(sc_counts))
  rm(sc_counts); gc()

  saveRDS(list(sc_dat = sc_dat, cell_type = cell_type, cell_state = cell_state),
          ref_rds)
}

## ============================================================================
## 2. 读取 TCGA star_counts bulk
## ============================================================================
if (file.exists(bulk_rds)) {
  message("读取已保存的 bulk 矩阵...")
  bk_dat <- readRDS(bulk_rds)
} else {
  f_counts <- list.files(TCGA_DIR, pattern = "star_counts.*\\.tsv\\.gz$",
                         full.names = TRUE, recursive = TRUE)
  f_probe  <- list.files(TCGA_DIR, pattern = "probemap", full.names = TRUE,
                         recursive = TRUE)
  stopifnot(length(f_counts) >= 1, length(f_probe) >= 1)
  dt  <- fread(f_counts[1])
  ens <- sub("\\..*$", "", dt[[1]])
  pm  <- fread(f_probe[1])
  sym <- pm$gene[match(ens, sub("\\..*$", "", pm$id))]
  keep <- !is.na(sym) & sym != ""
  expr <- as.matrix(dt[keep, -1, with = FALSE])
  rownames(expr) <- sym[keep]
  expr <- rowsum(expr, rownames(expr))
  bk_dat <- t(expr)
  saveRDS(bk_dat, bulk_rds)
  rm(dt, expr); gc()
}
message("bulk 矩阵: ", nrow(bk_dat), " 样本 x ", ncol(bk_dat), " 基因")

## ============================================================================
## 3. 运行 BayesPrism（跳过离群过滤版）
## ============================================================================
prism_rds <- file.path(OUT_DIR,
                       "08_prism_result.rds")

if (file.exists(prism_rds)) {
  bp <- readRDS(prism_rds)
} else {
  ## 状态→大类：每个亚型强制归到众数大类（保证树状结构）
  type_by_state <- sapply(split(cell_type, cell_state),
                          function(v) names(sort(table(v), decreasing = TRUE))[1])
  cell_type <- unname(type_by_state[cell_state])
  names(cell_type) <- names(cell_state)
  message("大类: ", paste(unique(cell_type), collapse = ", "))

  ## 直接取共同基因
  common_g <- intersect(colnames(sc_dat), colnames(bk_dat))
  message("去卷积使用的共同基因数: ", length(common_g))
  sc_f <- sc_dat[, common_g, drop = FALSE]
  bk_f <- bk_dat[, common_g, drop = FALSE]

  prism <- new.prism(reference = sc_f, mixture = bk_f,
                     input.type = "count.matrix",
                     cell.type.labels = cell_type,
                     cell.state.labels = cell_state,
                     key = "Tumor",
                     outlier.cut = 0.01, outlier.fraction = 0.1)
  message("开始 Gibbs 采样（1-3 小时），期间不要关闭 R ...")
  bp <- run.prism(prism = prism, n.core = 4)
  saveRDS(bp, prism_rds)
}
## ============================================================================
## 4. 提取 theta 并合并生存数据
## ============================================================================
theta <- as.data.frame(get.fraction(bp = bp, which.theta = "final",
                                    state.or.type = "state"))
message("theta: ", nrow(theta), " 样本 x ", ncol(theta), " 细胞状态")

f_surv <- list.files(TCGA_DIR, pattern = "survival.*\\.tsv(\\.gz)?$",
                     full.names = TRUE, recursive = TRUE)
surv <- as.data.frame(fread(f_surv[1]))
sid_col <- intersect(c("sample","_sample_id","submitter_id.samples"),
                     colnames(surv))[1]
rownames(surv) <- surv[[sid_col]]

common_s <- intersect(rownames(theta), rownames(surv))
theta <- theta[common_s, , drop = FALSE]
surv  <- surv[common_s, , drop = FALSE]
dat <- as.data.frame(cbind(theta, surv[, c("OS","OS.time")]))
dat$OS.time <- dat$OS.time / 30.44
dat <- dat[!is.na(dat$OS) & !is.na(dat$OS.time) & dat$OS.time > 0, ]

f_cli <- list.files(TCGA_DIR, pattern = "clinical.*\\.tsv(\\.gz)?$",
                    full.names = TRUE, recursive = TRUE)
cli <- as.data.frame(fread(f_cli[1]))
cid <- intersect(c("sample","_sample_id","submitter_id.samples"),
                 colnames(cli))[1]
dg_col <- intersect(c("primary_diagnosis","histological_type"), colnames(cli))[1]
dg  <- cli[[dg_col]][match(rownames(dat), cli[[cid]])]
dat$cohort <- ifelse(grepl("Glioblastoma", dg, ignore.case = TRUE), "GBM", "LGG")
message("最终分析样本: ", nrow(dat), " (GBM=", sum(dat$cohort=="GBM", na.rm=TRUE),
        ", LGG=", sum(dat$cohort=="LGG", na.rm=TRUE), ")")

## ============================================================================
## 5. theta × OS 单因素 Cox
## ============================================================================
subtypes <- colnames(theta)

run_cox <- function(df, label) {
  res <- lapply(subtypes, function(s) {
    x <- df[[s]]
    if (sd(x, na.rm = TRUE) == 0) return(NULL)
    fit <- tryCatch(coxph(Surv(OS.time, OS) ~ x, data = df),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    sfit <- summary(fit)
    data.frame(subtype = s, HR = sfit$conf.int[1, "exp(coef)"],
               lower = sfit$conf.int[1, "lower .95"],
               upper = sfit$conf.int[1, "upper .95"],
               p = sfit$coefficients[1, "Pr(>|z|)"], analysis = label)
  })
  rbindlist(res[!sapply(res, is.null)])
}

res_gbm <- run_cox(dat[dat$cohort == "GBM", ], "GBM_only")

res_all <- lapply(subtypes, function(s) {
  x <- dat[[s]]
  if (sd(x, na.rm = TRUE) == 0) return(NULL)
  fit <- tryCatch(coxph(Surv(OS.time, OS) ~ x + cohort, data = dat),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  sfit <- summary(fit)
  data.frame(subtype = s, HR = sfit$conf.int[1, "exp(coef)"],
             lower = sfit$conf.int[1, "lower .95"],
             upper = sfit$conf.int[1, "upper .95"],
             p = sfit$coefficients[1, "Pr(>|z|)"], analysis = "All_adj")
})
res_all <- rbindlist(res_all[!sapply(res_all, is.null)])

res <- rbind(res_gbm, res_all)
res$FDR <- p.adjust(res$p, method = "BH")
res <- res[order(res$analysis, res$p), ]
fwrite(res, file.path(OUT_DIR, "08_prism_theta_cox.csv"))

cat("\n========== FDR < 0.05 ==========\n")
print(as.data.frame(res[res$FDR < 0.05, ]))

message("08 完成。结果目录: ", OUT_DIR)
