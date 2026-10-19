#!/usr/bin/env Rscript
# ============================================================================
# Script      : 30_TF_activity_LDAM.R
# Purpose     : Transcription-factor activity inference in LDAM-high vs LDAM-low microglia (decoupleR + DoRothEA, ULM scores).
# Manuscript  : Supplementary analyses (upstream regulators)
# Input       : results/06_GSE190452_annotated_full.rds
# Output      : results/30_* TF activity tables/figures
# Dependencies: Seurat, UCell, decoupleR, dorothea, ggplot2, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================
# ===========================================================================
# 30_TF_activity_LDAM.R
# 收官分析 5：LDAM-high 小胶质细胞的转录因子活性推断
# 问题：谁在驱动 LDAM 程序？（转录层的上游调控者）
# 方法：decoupleR + DoRothEA 调控子库 + ULM 活性打分
#       （SCENIC 的轻量等效方案；调控子库经人工精选，A-C 置信度）
# 数据：GSE190452 小胶质细胞（results/06_GSE190452_annotated_full.rds）
# 分组：LDAM-high vs LDAM-low（与脚本 28 相同口径）
# 依赖：Seurat, decoupleR, dorothea, tidyverse, ggplot2
#   install.packages("decoupleR"); BiocManager::install("dorothea")
# ===========================================================================
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# install.packages("decoupleR"); BiocManager::install("dorothea")
suppressPackageStartupMessages({
  library(Seurat)
  library(decoupleR)
  library(dorothea)
  library(tidyverse)
  library(ggplot2)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/30_TF_activity", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. 读小胶质细胞对象
# ---------------------------------------------------------------------------
obj <- readRDS("results/06_GSE190452_annotated_full.rds")
ann_cols <- c("celltype", "cell_type", "annotation", "cluster")
ann_col <- intersect(ann_cols, colnames(obj@meta.data))[1]
mg <- subset(obj, cells = rownames(obj@meta.data)[grepl("micro|MG", obj@meta.data[[ann_col]], ignore.case = TRUE)])
message("小胶质细胞数: ", ncol(mg))

# ---------------------------------------------------------------------------
# 2. LDAM 分层（与脚本 28 相同：LDAM_ref_UCell > 中位数）
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(UCell))
LDAM_ref <- c("PLIN2","DGAT1","DGAT2","ACAT1","ACSL1","CD36","FABP3","HILPDA","SOAT1","LPIN1")
mg <- AddModuleScore_UCell(mg, features = list(LDAM_ref = LDAM_ref), ncores = 1)
mg$LDAM_group <- ifelse(mg$LDAM_ref_UCell > median(mg$LDAM_ref_UCell), "LDAM_high", "LDAM_low")
print(table(mg$LDAM_group))

# ---------------------------------------------------------------------------
# 3. 伪bulk 表达（每样本 × LDAM 组），避免单细胞零膨胀伪相关
# ---------------------------------------------------------------------------
mg <- JoinLayers(mg)

expr <- as.matrix(GetAssayData(mg, assay = "RNA", layer = "data"))
cell_sample <- mg@meta.data
sid_col <- intersect(c("orig.ident", "sample", "patient", "Sample"), colnames(cell_sample))[1]
message("样本列: ", sid_col)
pb_key <- paste(cell_sample[[sid_col]], mg$LDAM_group, sep = "__")
pb <- rowsum(t(expr), group = pb_key)          # 行 = 样本__组, 列 = 基因
pb <- t(pb)
message("伪bulk 维度: ", nrow(pb), " 基因 x ", ncol(pb), " 样本组")
# 过滤：总表达过低的伪bulk 列
pb <- pb[, colSums(pb) > 0]

# log-CPM 标准化（供 ULM 使用）
pb_norm <- t(t(pb) / colSums(pb) * 1e6)
pb_norm <- log2(pb_norm + 1)

# ---------------------------------------------------------------------------
# 4. DoRothEA 调控子 + decoupleR ULM 活性
# ---------------------------------------------------------------------------
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C")) %>%
  select(source = tf, target, mor)

acts <- run_ulm(pb_norm, regulons, .source = "source", .target = "target", .mor = "mor")
write_csv(acts, "results/30_TF_activity/30_TF_activity_all_pseudobulk.csv")

# ---------------------------------------------------------------------------
# 5. LDAM-high vs low 差异活性（按样本配对无法保证，用非配对 Wilcoxon）
# ---------------------------------------------------------------------------
acts$group <- sub(".*__", "", acts$condition)
acts$sample <- sub("__.*", "", acts$condition)

tf_res <- acts %>% group_by(source) %>%
  summarise(delta = mean(score[group == "LDAM_high"], na.rm = TRUE) -
                    mean(score[group == "LDAM_low"], na.rm = TRUE),
            P = tryCatch(wilcox.test(score[group == "LDAM_high"],
                                     score[group == "LDAM_low"])$p.value,
                         error = function(e) NA_real_),
            .groups = "drop")
tf_res$FDR <- p.adjust(tf_res$P, "BH")
tf_res <- tf_res %>% arrange(P)
print(tf_res, n = 30)
write_csv(tf_res, "results/30_TF_activity/30_TF_differential_activity.csv")

# ---------------------------------------------------------------------------
# 6. 图：top 25 差异 TF
# ---------------------------------------------------------------------------
top25 <- tf_res %>% filter(!is.na(P)) %>% slice_head(n = 25) %>%
  mutate(source = fct_reorder(source, delta))
p <- ggplot(top25, aes(delta, source, fill = FDR < 0.05)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#C9828B", "FALSE" = "#B9C4D4")) +
  labs(title = "Differential TF activity: LDAM-high vs LDAM-low microglia (GSE190452)",
       subtitle = "DoRothEA regulons (A-C confidence), ULM scores on pseudobulk",
       x = "activity difference (high - low)", y = NULL, fill = "FDR<0.05")
ggsave("results/30_TF_activity/30_TF_activity_top25.pdf", p, width = 7, height = 7)

message("完成。判读：")
message("- 预期阳性：PPARG/SREBF1/SREBF2（脂质合成）、CEBPA/CEBPB（LDAM 经典驱动）")
message("- 铁死亡应激侧：NFE2L2(NRF2)/ATF4/HIF1A 若显著，补上转录层的应激证据")
message("- 若 PPARG 显著：与 32 的 PPAR 激动剂命中形成「机制-药物」呼应")
