# ============================================================================
# Script      : 08_blood_GSE143272_mouse_GSE88992.R
# Purpose     : Part A: GSE143272 peripheral blood (drug-naive epilepsy vs healthy) preprocessing for 04; Part B: mouse KA hippocampus array (GSE88992) cross-species validation.
# Manuscript  : Methods 2.5 / Results 3.7
# Input       : GSE143272_non-normalized.txt.gz; GSE88992_RAW.tar; genesets/
# Output      : normalized blood matrix for 04; mouse validation tables/figures
# Dependencies: GSVA, affy, biomaRt, limma, mouse4302.db
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================

# ==============================================================================
# 08_blood_GSE143272_mouse_GSE88992.R
# Part A：GSE143272 人外周血（未用药癫痫 34 vs 健康对照 50）——外周标志物验证队列
#         输出归一化矩阵供脚本 04 做外周血 ROC 验证
# Part B：GSE88992 小鼠 KA 模型海马芯片 —— 动物模型层面旁证（与阶段 3 呼应）
# ==============================================================================

# ============================ Part A：GSE143272 外周血（全程不联网） ============================

## A.1 读取非归一化矩阵
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root

f143 <- list.files(BASE_DIR, pattern = "GSE143272.*\\.txt\\.gz$", full.names = TRUE)[1]
message("读取: ", f143)
raw <- fread(f143)

## A.2 提取信号列与 Detection P 值列（该文件信号列即样本名）
det_cols <- grep("Detection", names(raw), value = TRUE)
sig_cols <- setdiff(names(raw), c(names(raw)[1], det_cols))
cat("信号列（样本）数：", length(sig_cols), "\n")        # 应为 142

expr <- as.matrix(raw[, ..sig_cols])
storage.mode(expr) <- "numeric"
rownames(expr) <- raw[[1]]
colnames(expr) <- iconv(colnames(expr), to = "UTF-8", sub = "")
det_mat <- as.matrix(raw[, ..det_cols])
storage.mode(det_mat) <- "numeric"

## A.3 neqc 背景校正 + 归一化 + 探针 → symbol（GPL10558 本地注释）
expr_neqc <- neqc(expr, detection.p = det_mat)

gpl_file <- file.path(BASE_DIR, "GSE143272", "GPL10558.annot.gz")
if (!file.exists(gpl_file)) {
  stop("请先用浏览器下载：\n",
       "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz\n",
       "保存到：", gpl_file, " 后重新运行。")
}
gpl_tab   <- Table(getGEO(filename = gpl_file))
sym_col   <- grep("symbol", names(gpl_tab), ignore.case = TRUE, value = TRUE)[1]
id_col_gpl <- grep("^ID$", names(gpl_tab), value = TRUE)[1]
probe2gene <- setNames(as.character(gpl_tab[[sym_col]]),
                       as.character(gpl_tab[[id_col_gpl]]))

sym_hit <- probe2gene[rownames(expr_neqc)]
keep <- !is.na(sym_hit) & sym_hit != ""
expr_neqc <- expr_neqc[keep, ]; sym_hit <- sym_hit[keep]
sym_hit <- str_split_fixed(sym_hit, " /// ", 2)[, 1]
cat("成功注释探针：", length(sym_hit), "\n")

expr_gene <- avereps(expr_neqc, ID = sym_hit)
cat("基因水平矩阵：", nrow(expr_gene), "基因 x", ncol(expr_gene), "样本\n")

## A.4 分组（useBytes=TRUE 绕开 "Drug-naïve" 特殊字符的编码问题）
nm <- colnames(expr_gene)
grp <- rep(NA_character_, length(nm))
grp[grepl("^Healthy", nm, useBytes = TRUE)] <- "Control"
grp[grepl("^Drug",    nm, useBytes = TRUE)] <- "Epilepsy"
keep_s <- !is.na(grp)
expr_use <- expr_gene[, keep_s]
meta_use <- data.frame(sample = colnames(expr_use), group = grp[keep_s])
rownames(meta_use) <- meta_use$sample
meta_use$group <- factor(meta_use$group, levels = c("Control", "Epilepsy"))
print(table(meta_use$group))          # ！！应为 Control≈50, Epilepsy≈34
## 清洗列名中的非法编码字节（"Drug-na\xefve" → "Drug-nave"）
colnames(expr_use) <- iconv(colnames(expr_use), to = "UTF-8", sub = "")
meta_use$sample <- colnames(expr_use)
rownames(meta_use) <- meta_use$sample

## 然后从 fwrite 那两行开始，把 A.4 剩余部分 + A.5 + A.6 重新跑一遍
stopifnot(all(c("Control", "Epilepsy") %in% levels(droplevels(meta_use$group))))

fwrite(data.frame(gene = rownames(expr_use), expr_use, check.names = FALSE),
       file.path(OUT_DIR, "08_GSE143272_normalized.csv"))
fwrite(meta_use, file.path(OUT_DIR, "08_GSE143272_sample_meta.csv"))

## A.5 外周血 LDAM / 铁死亡 ssGSEA 评分 + 组间比较
library(GSVA)
sc <- gsva(gsvaParam(as.matrix(expr_use),
                     list(LDAM = intersect(LDAM_genes, rownames(expr_use)),
                          Ferroptosis = intersect(Ferroptosis_genes, rownames(expr_use)))),
           verbose = FALSE)
fwrite(data.frame(sample = colnames(sc), t(sc)),
       file.path(OUT_DIR, "08_GSE143272_pathway_scores.csv"))

df_b <- data.frame(t(sc), group = meta_use$group)
p_b <- df_b |>
  tidyr::pivot_longer(c(LDAM, Ferroptosis), names_to = "pathway", values_to = "score") |>
  ggplot(aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~pathway, scales = "free") +
  scale_fill_manual(values = c("Control" = "#7BA7BC", "Epilepsy" = "#C0392B")) +
  stat_compare_means(method = "wilcox.test") +
  labs(title = "GSE143272 blood: LDAM / Ferroptosis ssGSEA score") + theme_bw()
ggsave(file.path(FIG_DIR, "08_blood_LDam_Ferroptosis_score.pdf"), p_b, width = 8, height = 5)

## 组间 Wilcoxon 检验直接打印（图上看不清时以这里为准）
wilcox.test(sc["LDAM", ] ~ meta_use$group) |> print()
wilcox.test(sc["Ferroptosis", ] ~ meta_use$group) |> print()

## A.6 差异分析（Epilepsy vs Control，limma）——供 04 号脚本交叉验证用
design <- model.matrix(~ meta_use$group)
fit <- eBayes(lmFit(expr_use, design))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
res_blood <- data.table(gene = rownames(tt), log2FoldChange = tt$logFC,
                        pvalue = tt$P.Value, padj = tt$adj.P.Val)
res_blood[, direction := fifelse(padj < 0.05 & log2FoldChange >  0.5, "Up",
                                 fifelse(padj < 0.05 & log2FoldChange < -0.5, "Down", "NS"))]
fwrite(res_blood, file.path(OUT_DIR, "08_GSE143272_DEG_Epilepsy_vs_Control.csv"))
cat("== GSE143272 血液 DEG 数 ==\n"); print(table(res_blood$direction))

cat("\n== Part A 完成！==\n")
# ============================ Part B：GSE88992 小鼠 KA 模型旁证 ============================
## B.1 读取 CEL 文件（Affymetrix Mouse 430 2.0 芯片）
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# BiocManager::install(c("affy", "mouse4302.db"), ask = FALSE, update = FALSE)
library(limma)
untar_dir <- file.path(BASE_DIR, "GSE88992")
if (!dir.exists(untar_dir)) {
  dir.create(untar_dir)
  untar(file.path(BASE_DIR, "GSE88992_RAW.tar"), exdir = untar_dir)
}
cel_files <- list.files(untar_dir, pattern = "\\.CEL(\\.gz)?$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
print(cel_files)   # 17 个混池样本：KA 8 + saline 9

if (length(cel_files) > 0 && requireNamespace("affy", quietly = TRUE)) {
  library(affy)
  ab <- ReadAffy(filenames = cel_files)
  eset_m <- rma(ab)
  expr_m <- exprs(eset_m)

  # 探针 → 小鼠 symbol（mouse4302.db）
  if (requireNamespace("mouse4302.db", quietly = TRUE)) {
    library(mouse4302.db)
    sym <- AnnotationDbi::mapIds(mouse4302.db, keys = rownames(expr_m),
                                 keytype = "PROBEID", column = "SYMBOL",
                                 multiVals = "first")
  } else {
    # 备用：GEO 平台注释
    gpl_m <- Table(getGEO("GPL1261", AnnotGPL = FALSE))
    scol <- grep("symbol", names(gpl_m), ignore.case = TRUE, value = TRUE)[1]
    icol <- grep("^ID$", names(gpl_m), value = TRUE)[1]
    sym <- setNames(gpl_m[[scol]], gpl_m[[icol]])[rownames(expr_m)]
  }
  keep <- !is.na(sym) & sym != ""
  expr_mg <- aggregate(expr_m[keep, ], by = list(sym[keep]), FUN = mean)
  rn <- expr_mg[[1]]; expr_mg <- as.matrix(expr_mg[, -1]); rownames(expr_mg) <- rn

  # 分组：文件名含 KA/kainate → KA，saline/control → Saline（按打印文件名核对）
  group_m <- ifelse(str_detect(tolower(colnames(expr_mg)), "ka|kainate"), "KA", "Saline")
  meta_m <- data.frame(sample = colnames(expr_mg), group = group_m)
  print(table(meta_m$group))

  # 小鼠基因 → 人同源基因（biomaRt）；失败则直接用小鼠基因做描述性分析
  ldam_m <- LDAM_genes; ferr_m <- Ferroptosis_genes
  tryCatch({
    library(biomaRt)
    hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    conv <- getLDS(attributes = "hgnc_symbol", filters = "hgnc_symbol",
                   values = unique(c(LDAM_genes, Ferroptosis_genes)), mart = hs,
                   attributesL = "mgi_symbol", martL = mm)
    ldam_m <- unique(conv$MGI.symbol[conv$HGNC.symbol %in% LDAM_genes])
    ferr_m <- unique(conv$MGI.symbol[conv$HGNC.symbol %in% Ferroptosis_genes])
    message("同源转换完成: LDAM ", length(ldam_m), " / 铁死亡 ", length(ferr_m), " 个小鼠基因")
  }, error = function(e) message("biomaRt 失败，跳过同源转换: ", conditionMessage(e)))

  ## B.2 差异分析（limma）+ 小鼠 LDAM/铁死亡评分
  design <- model.matrix(~ 0 + factor(meta_m$group))
  colnames(design) <- c("KA", "Saline")
  fit <- lmFit(expr_mg, design) |>
    contrasts.fit(makeContrasts(KA - Saline, levels = design)) |> eBayes()
  de_m <- topTable(fit, number = Inf, adjust = "BH") |>
    tibble::rownames_to_column("gene")
  fwrite(de_m, file.path(OUT_DIR, "08_GSE88992_mouse_DEG_KA_vs_Saline.csv"))

  library(GSVA)
  sc_m <- gsva(gsvaParam(expr_mg, list(LDAM = ldam_m, Ferroptosis = ferr_m)))
  df_m <- data.frame(t(sc_m), group = meta_m$group)
  p_m <- df_m |>
    tidyr::pivot_longer(c(LDAM, Ferroptosis), names_to = "pathway", values_to = "score") |>
    ggplot(aes(group, score, fill = group)) +
    geom_boxplot() + geom_point() + facet_wrap(~pathway, scales = "free") +
    scale_fill_manual(values = c("Saline" = "#7BA7BC", "KA" = "#C0392B")) +
    stat_compare_means(method = "wilcox.test") +
    labs(title = "GSE88992 小鼠 KA 模型海马：LDAM / 铁死亡评分") + theme_bw()
  ggsave(file.path(FIG_DIR, "08_mouse_LDam_Ferroptosis_score.pdf"), p_m,
         width = 8, height = 5)
} else {
  message("未找到 CEL 文件或 affy 包未安装，跳过小鼠分析（BiocManager::install('affy')）")
}

message("=== 08_blood_GSE143272_mouse_GSE88992.R 运行完毕 ===")
## 人类基因名 → 小鼠写法（首字母大写其余小写）
to_mouse <- function(g) paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
ldam_m <- intersect(to_mouse(LDAM_genes), rownames(expr_mg))
ferr_m <- intersect(to_mouse(Ferroptosis_genes), rownames(expr_mg))
cat("匹配到小鼠基因：LDAM", length(ldam_m), "个；铁死亡", length(ferr_m), "个\n")

## 重新打分 + 出图
library(GSVA)
sc_m <- gsva(gsvaParam(expr_mg, list(LDAM = ldam_m, Ferroptosis = ferr_m)))
fwrite(data.frame(sample = colnames(sc_m), t(sc_m), group = meta_m$group),
       file.path(OUT_DIR, "08_GSE88992_pathway_scores.csv"))

df_m <- data.frame(t(sc_m), group = meta_m$group)
p_m <- df_m |>
  tidyr::pivot_longer(c(LDAM, Ferroptosis), names_to = "pathway", values_to = "score") |>
  ggplot(aes(group, score, fill = group)) +
  geom_boxplot() + geom_point() + facet_wrap(~pathway, scales = "free") +
  scale_fill_manual(values = c("Saline" = "#7BA7BC", "KA" = "#C0392B")) +
  stat_compare_means(method = "wilcox.test") +
  labs(title = "GSE88992 mouse KA hippocampus: LDAM / Ferroptosis") + theme_bw()
ggsave(file.path(FIG_DIR, "08_mouse_LDam_Ferroptosis_score.pdf"), p_m, width = 8, height = 5)

## 控制台直接看检验结果
wilcox.test(sc_m["LDAM", ] ~ meta_m$group) |> print()
wilcox.test(sc_m["Ferroptosis", ] ~ meta_m$group) |> print()
