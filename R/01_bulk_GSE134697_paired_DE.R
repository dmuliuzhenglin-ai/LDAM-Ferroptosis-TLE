# ============================================================================
# Script      : 01_bulk_GSE134697_paired_DE.R
# Purpose     : Paired DESeq2 differential expression in GSE134697 (17 mTLE hippocampus vs autologous neocortex).
# Manuscript  : Methods 2.2 / Results 3.1
# Input       : GSE134697_CountMatrix.tsv.gz
# Output      : results/01_* DE tables; figures (volcano etc.)
# Dependencies: DESeq2, data.table, dplyr, ggplot2, ggrepel, org.Hs.eg.db, stringr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

# ==============================================================================
# 01_bulk_GSE134697_paired_DE.R
# GSE134697：17 例 mTLE 患者配对（海马 vs 自体新皮层）+ 2 例健康新皮层对照
# 设计：DESeq2 配对设计（~ patient + tissue），海马 vs 皮层
# 数据文件：GSE134697_CountMatrix.tsv.gz
# ==============================================================================

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results"); FIG_DIR <- file.path(BASE_DIR, "figures")
dir.create(OUT_DIR, showWarnings = FALSE); dir.create(FIG_DIR, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table); library(DESeq2); library(dplyr); library(ggplot2)
  library(ggrepel); library(stringr)
})

## ---------- 1.1 读取 counts 矩阵 ----------
# GSE134697_CountMatrix.tsv.gz 为 featureCounts 原始输出：
# 第 1 列 Geneid，第 2-6 列为 Chr/Start/End/Strand/Length 注释，其后为样本计数列
# 样本列命名：CHKJ-H1~H17=患者海马，CHKJ-C1~C17=患者自体皮层（与 H 编号一一配对），
#             CTRL-00001/00002=健康对照皮层
cnt <- fread(file.path(BASE_DIR, "GSE134697_CountMatrix.tsv.gz"))
message("矩阵维度: ", nrow(cnt), " 行 x ", ncol(cnt), " 列")

sample_cols <- grep("\\.bam$", names(cnt), value = TRUE)
expr <- as.data.frame(cnt[, ..sample_cols])
rownames(expr) <- cnt[[1]]                    # 第 1 列 Geneid

# 简化列名为 CHKJ-H1 / CHKJ-C1 / CTRL-00001 短名
short_names <- str_extract(basename(sample_cols), "CHKJ-[CH]\\d+|CTRL-\\d+")
colnames(expr) <- short_names

expr[] <- lapply(expr, as.numeric)
expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", ]
expr <- round(rowsum(as.matrix(expr), group = rownames(expr), reorder = FALSE))  # 重复基因合并
storage.mode(expr) <- "integer"

# 若基因 ID 为 ENSEMBL（ENSG 开头）则转为 symbol
if (mean(grepl("^ENSG", rownames(expr))) > 0.5) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = sub("\\..*$", "", rownames(expr)),
                               keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  keep <- !is.na(sym)
  expr <- rowsum(expr[keep, ], group = sym[keep], reorder = FALSE)
}

## ---------- 1.2 样本分组 ----------
meta <- data.frame(
  sample  = colnames(expr),
  tissue  = case_when(
    str_detect(colnames(expr), "^CHKJ-H") ~ "Hippocampus",
    str_detect(colnames(expr), "^CHKJ-C") ~ "Neocortex",
    str_detect(colnames(expr), "^CTRL")   ~ "Control_neocortex"),
  patient = ifelse(str_detect(colnames(expr), "^CHKJ"),
                   paste0("P", str_extract(colnames(expr), "\\d+")), NA)
)
print(table(meta$tissue))
print(meta)
# ！！请核对上表：应得到 Hippocampus=17, Neocortex=17, Control_neocortex=2
# 且 P1-P17 每个患者都有 H 和 C 两行（配对关系）

## ---------- 1.3 配对差异分析：患者海马 vs 患者自体皮层 ----------
pair_idx <- meta$tissue %in% c("Hippocampus", "Neocortex")
pair_meta <- meta[pair_idx, ]
pair_meta <- pair_meta[!is.na(pair_meta$patient) &
                         pair_meta$patient %in% names(which(table(pair_meta$patient) == 2)), ]
pair_meta$tissue  <- factor(pair_meta$tissue, levels = c("Neocortex", "Hippocampus"))
pair_meta$patient <- factor(pair_meta$patient)

dds <- DESeqDataSetFromMatrix(countData = expr[, pair_meta$sample],
                              colData   = pair_meta,
                              design    = ~ patient + tissue)   # 配对设计
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]                  # 低表达过滤
dds <- DESeq(dds)
res <- results(dds, contrast = c("tissue", "Hippocampus", "Neocortex"))

de_df <- as.data.frame(res) |>
  tibble::rownames_to_column("gene") |>
  arrange(padj) |>
  mutate(group = case_when(padj < 0.05 & log2FoldChange >  1 ~ "Up",
                           padj < 0.05 & log2FoldChange < -1 ~ "Down",
                           TRUE ~ "NS"))
fwrite(de_df, file.path(OUT_DIR, "01_GSE134697_DEG_Hippocampus_vs_Neocortex.csv"))
message("显著上调: ", sum(de_df$group == "Up", na.rm = TRUE),
        " | 显著下调: ", sum(de_df$group == "Down", na.rm = TRUE))

## ---------- 1.4 火山图 ----------
p_vol <- ggplot(de_df, aes(log2FoldChange, -log10(padj), color = group)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up" = "#C0392B", "Down" = "#2C5F8A", "NS" = "grey80")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
  labs(title = "GSE134697 配对分析：TLE 海马 vs 自体新皮层",
       x = "log2 Fold Change", y = "-log10(padj)") +
  theme_bw()
ggsave(file.path(FIG_DIR, "01_volcano_GSE134697.pdf"), p_vol, width = 7, height = 6)

## ---------- 1.5 归一化矩阵导出（供 GSVA / 后续脚本使用） ----------
vsd <- vst(dds, blind = FALSE)
norm_mat <- assay(vsd)
fwrite(data.frame(gene = rownames(norm_mat), norm_mat),
       file.path(OUT_DIR, "01_GSE134697_vst_normalized.csv"))
fwrite(pair_meta, file.path(OUT_DIR, "01_GSE134697_sample_meta.csv"))

message("=== 01_bulk_GSE134697_paired_DE.R 运行完毕 ===")
