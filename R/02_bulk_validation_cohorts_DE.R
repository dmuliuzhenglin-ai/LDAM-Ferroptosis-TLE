# ============================================================================
# Script      : 02_bulk_validation_cohorts_DE.R
# Purpose     : Differential expression in two validation cohorts: GSE190451 (human TLE cortex bulk) and GSE168375 (sorted astrocytes).
# Manuscript  : Methods 2.2 / Results 3.1
# Input       : GSE190451_TLE_samples.txt.gz; GSE168375_RAW.tar
# Output      : results/02_* DE tables and concordance checks
# Dependencies: DESeq2, GEOquery, data.table, dplyr, ggplot2, ggrepel, limma, org.Hs.eg.db, stringr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

# ==============================================================================
# 02_bulk_validation_cohorts_DE.R  【完整重写版】
# 两个带对照的验证队列：
#   Part A：GSE190451 —— 人 TLE 颞叶皮层 bulk RNA-seq，3 TLE vs 3 对照
#           （组织级验证队列，GSE190452 单细胞的配套 bulk）
#   Part B：GSE168375 —— 皮层纯化星形胶质细胞 RNA-seq（Control/Epilepsy，
#           Peritumor 瘤周样本剔除），星形胶质细胞层面验证队列
# 数据文件：GSE190451_TLE_samples.txt.gz（合并矩阵）、GSE168375_RAW.tar（逐样本）
# ==============================================================================
# ===== 先补这一步：读取 GSE190451 文件 =====
# [EN] Patch step: read the GSE190451 file first
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root

f451 <- list.files(BASE_DIR, pattern = "GSE190451.*\\.txt\\.gz$", full.names = TRUE)[1]
d451 <- fread(f451)
# [debug] print(names(d451))        # 应显示 Gene_id | Normal_1~3 | TLE_1~3 | Gene_name | Gene_description

# ===== 然后接着运行修正代码 =====
# [EN] Then run the corrected code below
sample_cols <- grep("^(Normal|TLE)_", names(d451), value = TRUE)
print(sample_cols)

exprA <- as.matrix(d451[, ..sample_cols])
storage.mode(exprA) <- "double"
sym <- d451$Gene_name
keep <- !is.na(sym) & sym != ""
exprA <- rowsum(exprA[keep, ], group = sym[keep], reorder = FALSE)
message("基因数: ", nrow(exprA), " | 样本数: ", ncol(exprA))
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results"); FIG_DIR <- file.path(BASE_DIR, "figures")

suppressPackageStartupMessages({
  library(data.table); library(DESeq2); library(limma); library(dplyr)
  library(ggplot2); library(ggrepel); library(stringr); library(GEOquery)
})

# 通用工具：基因 ID 自动识别并转 symbol
to_symbol <- function(mat) {
  ids <- rownames(mat)
  if (mean(grepl("^[0-9]+$", head(ids, 500))) > 0.5) {
    message("检测到 Entrez ID，转换为 symbol")
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids, keytype = "ENTREZID",
                                 column = "SYMBOL", multiVals = "first")
    keep <- !is.na(sym) & sym != ""
    mat <- rowsum(mat[keep, , drop = FALSE], group = sym[keep], reorder = FALSE)
  } else if (mean(grepl("^ENSG", head(ids, 500))) > 0.5) {
    message("检测到 ENSEMBL ID，转换为 symbol")
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = sub("\\..*$", "", ids),
                                 keytype = "ENSEMBL", column = "SYMBOL",
                                 multiVals = "first")
    keep <- !is.na(sym) & sym != ""
    mat <- rowsum(mat[keep, , drop = FALSE], group = sym[keep], reorder = FALSE)
  }
  mat
}

# 通用工具：火山图
plot_volcano <- function(de_df, title, outfile) {
  p <- ggplot(de_df, aes(log2FoldChange, -log10(padj), color = group)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "#C0392B", "Down" = "#2C5F8A", "NS" = "grey80")) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(padj)") + theme_bw()
  ggsave(outfile, p, width = 7, height = 6)
}

# 通用工具：counts 用 DESeq2 / 非整数用 limma，统一输出标准 DEG 表
run_de <- function(expr, meta, contrast_levels, tag) {
  is_counts <- all(expr == floor(expr), na.rm = TRUE)
  if (is_counts) {
    message(tag, ": 检测到整数 counts，使用 DESeq2")
    dds <- DESeqDataSetFromMatrix(round(expr), meta, design = ~ group)
    dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("group", contrast_levels[1], contrast_levels[2]))
    de_df <- as.data.frame(res) |> tibble::rownames_to_column("gene")
    norm_mat <- assay(vst(dds, blind = FALSE))
  } else {
    message(tag, ": 检测到非整数表达值，使用 limma（log2 转换）")
    log_expr <- log2(expr + 1)
    design <- model.matrix(~ meta$group)
    fit <- lmFit(log_expr, design) |> eBayes()
    de_df <- topTable(fit, coef = 2, number = Inf, adjust = "BH") |>
      tibble::rownames_to_column("gene") |>
      dplyr::rename(log2FoldChange = logFC, padj = adj.P.Val)
    norm_mat <- log_expr
  }
  de_df <- de_df |>
    arrange(padj) |>
    mutate(group = case_when(padj < 0.05 & log2FoldChange >  1 ~ "Up",
                             padj < 0.05 & log2FoldChange < -1 ~ "Down",
                             TRUE ~ "NS"))
  list(de = de_df, norm = norm_mat)
}

## ---------- Part A 读取（修正版） ----------
# 列结构：Gene_id | Normal_1~3 | TLE_1~3 | Gene_name | Gene_description
# 只取 6 个样本列，基因名直接用 Gene_name 列
sample_cols <- grep("^(Normal|TLE)_", names(d451), value = TRUE)
print(sample_cols)   # 应为 Normal_1/2/3 + TLE_1/2/3

exprA <- as.matrix(d451[, ..sample_cols])
storage.mode(exprA) <- "double"

sym <- d451$Gene_name
keep <- !is.na(sym) & sym != ""
exprA <- rowsum(exprA[keep, ], group = sym[keep], reorder = FALSE)  # 重复基因合并
message("基因数: ", nrow(exprA), " | 样本数: ", ncol(exprA))

# ==================== Part B：GSE168375（星形胶质细胞层面验证队列） ====================
message("========== Part B: GSE168375 ==========")

untar_dir <- file.path(BASE_DIR, "GSE168375")
if (!dir.exists(untar_dir)) {
  dir.create(untar_dir)
  untar(file.path(BASE_DIR, "GSE168375_RAW.tar"), exdir = untar_dir)
}
files <- list.files(untar_dir, pattern = "\\.txt\\.gz$", recursive = TRUE,
                    full.names = TRUE)
print(basename(files))       # 应约 49 个 GSM 开头的 counts 文件

read_one <- function(f) {
  d <- fread(f)
  v <- as.numeric(d[[ncol(d)]])            # 最后一列 = counts
  names(v) <- as.character(d[[1]])         # 第一列 = 基因 ID
  tapply(v, names(v), sum)                 # 同文件内基因去重
}
count_list <- lapply(files, read_one)
count_list <- count_list[sapply(count_list, length) > 1000]
message("有效样本文件数: ", length(count_list))

all_genes <- unique(unlist(lapply(count_list, names)))
exprB <- do.call(cbind, lapply(count_list, function(v) v[all_genes]))
rownames(exprB) <- all_genes
colnames(exprB) <- str_extract(basename(files), "GSM\\d+")
exprB[is.na(exprB)] <- 0
exprB <- to_symbol(exprB)
message("基因数: ", nrow(exprB), " | 样本数: ", ncol(exprB))

# 分组：GEO 标题前缀 Control_/Epilepsy_，剔除 Peritumor
gse <- getGEO("GSE168375", GSEMatrix = TRUE, getGPL = FALSE)
pheno <- pData(gse[[1]])
titles <- setNames(pheno$title, pheno$geo_accession)[colnames(exprB)]
print(table(sub("_.*$", "", titles)))      # 应显示 Control / Epilepsy / Peritumor

groupB <- case_when(
  str_detect(titles, "^Control")  ~ "Control",
  str_detect(titles, "^Epilepsy") ~ "Epilepsy",
  TRUE                            ~ NA_character_)
metaB <- data.frame(sample = colnames(exprB), group = groupB) |> filter(!is.na(group))
metaB$group <- factor(metaB$group, levels = c("Control", "Epilepsy"))
exprB <- exprB[, metaB$sample]
print(table(metaB$group))    # ！！应约为 Control=12, Epilepsy=31

resB <- run_de(exprB, metaB, c("Epilepsy", "Control"), "GSE168375")
fwrite(resB$de, file.path(OUT_DIR, "02_GSE168375_DEG_Epilepsy_vs_Control.csv"))
fwrite(data.frame(gene = rownames(resB$norm), resB$norm),
       file.path(OUT_DIR, "02_GSE168375_vst_normalized.csv"))
fwrite(metaB, file.path(OUT_DIR, "02_GSE168375_sample_meta.csv"))
plot_volcano(resB$de, "GSE168375：癫痫 vs 对照（星形胶质细胞）",
             file.path(FIG_DIR, "02_volcano_GSE168375.pdf"))
message("GSE168375 上调: ", sum(resB$de$group == "Up", na.rm = TRUE),
        " | 下调: ", sum(resB$de$group == "Down", na.rm = TRUE))
list.files(file.path(BASE_DIR, "results"))
## ========== 核对已写出的 02 结果文件 ==========
## [EN] Verify the 02 result files written to disk
library(data.table)

res_dir <- file.path(BASE_DIR, "results")

## ---- Part B: GSE168375 ----
resB <- fread(file.path(res_dir, "02_GSE168375_DEG_Epilepsy_vs_Control.csv"))
cat("== GSE168375 分组 ==\n"); print(table(fread(file.path(res_dir, "02_GSE168375_sample_meta.csv"))$group))
cat("== GSE168375 DEG 数 ==\n"); print(table(resB$direction))
cat("== GSE168375 前 10 个显著基因 ==\n")
# [debug] print(head(resB[order(padj), .(gene, log2FoldChange, padj, direction)], 10))

## ---- Part A: GSE190451 文件是否存在于别处 ----
cat("\n== 搜索 GSE190451 相关结果文件 ==\n")
print(list.files(res_dir, pattern = "190451", full.names = TRUE))
## ---------- 两队列差异方向一致性检查（重叠基因 log2FC 相关性） ----------
common <- intersect(resA$de$gene, resB$de$gene)
df_c <- data.frame(gene = common,
                   log2FC_GSE190451 = resA$de$log2FoldChange[match(common, resA$de$gene)],
                   log2FC_GSE168375 = resB$de$log2FoldChange[match(common, resB$de$gene)])
cor_c <- cor.test(df_c$log2FC_GSE190451, df_c$log2FC_GSE168375, method = "spearman")
message("两队列 log2FC 一致性: Spearman r = ", round(cor_c$estimate, 3),
        ", P = ", signif(cor_c$p.value, 3))

message("=== 02_bulk_validation_cohorts_DE.R 运行完毕 ===")

# 如果环境还在，也可以直接在控制台运行这几行快速查看：
table(metaA$group)        # Part A 分组：应 Control=3, TLE=3
table(metaB$group)        # Part B 分组：应约 Control=12, Epilepsy=31
sum(resA$de$group=="Up", na.rm=TRUE); sum(resA$de$group=="Down", na.rm=TRUE)   # GSE190451 上调/下调数
sum(resB$de$group=="Up", na.rm=TRUE); sum(resB$de$group=="Down", na.rm=TRUE)   # GSE168375 上调/下调数
