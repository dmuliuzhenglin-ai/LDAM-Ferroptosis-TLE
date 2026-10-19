#!/usr/bin/env Rscript
# ============================================================================
# Script      : 31_pseudotime_LDAM.R
# Purpose     : Monocle3 pseudotime trajectory from LDAM-low to LDAM-high microglia; dynamics of ferroptosis/exosome scores along the trajectory.
# Manuscript  : Supplementary analyses (trajectory)
# Input       : results/06_GSE190452_annotated_full.rds
# Output      : results/31_* trajectory tables/figures
# Dependencies: Seurat, UCell, ggplot2, monocle3, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================
# ===========================================================================
# 31_pseudotime_LDAM.R
# 收官分析 4：小胶质细胞 LDAM-low -> LDAM-high 拟时序轨迹
# 问题：LDAM 高低亚群是静态两群，还是同一病程轴上的连续演化？
#       铁死亡评分、外泌体机器评分沿轨迹如何变化？
# 数据：GSE190452 小胶质细胞（results/06_GSE190452_annotated_full.rds）
# 方法：Monocle3 拟时序 + 沿轨迹的评分动态图
# 依赖：Seurat, monocle3, UCell, tidyverse, ggplot2
# 注意：monocle3 在 Windows 安装若失败：
#   BiocManager::install(c("monocle3"))  或  remotes::install_github("cole-trapnell-lab/monocle3")
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(UCell)
  library(tidyverse)
  library(ggplot2)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/31_pseudotime", showWarnings = FALSE, recursive = TRUE)

gene_sets <- list(
  Exosome_uptake  = c("CD9","CD81","ICAM1","ITGB1","ITGAV","SDC1","SDC4","HSPG2"),
  Cargo_loading   = c("YBX1","SYNCRIP","HNRNPA2B1","FMR1","LIN28B","IGF2BP1"),
  Rab_secretory   = c("RAB27A","RAB27B","RAB11A","RAB11B","RAB35","RAB8A","RAB8B","RAB7A"),
  LDAM_ref        = c("PLIN2","DGAT1","DGAT2","ACAT1","ACSL1","CD36","FABP3","HILPDA","SOAT1","LPIN1"),
  Ferroptosis_ref = c("GPX4","ACSL4","LPCAT3","ALOX15","AIFM2","TFRC","FTH1","FTL","SLC7A11","NCOA4")
)

# ---------------------------------------------------------------------------
# 1. 读小胶质细胞对象 + UCell 打分
# ---------------------------------------------------------------------------
obj <- readRDS("results/06_GSE190452_annotated_full.rds")
ann_cols <- c("celltype", "cell_type", "annotation", "cluster")
ann_col <- intersect(ann_cols, colnames(obj@meta.data))[1]
mg <- subset(obj, cells = rownames(obj@meta.data)[grepl("micro|MG", obj@meta.data[[ann_col]], ignore.case = TRUE)])
message("小胶质细胞数: ", ncol(mg))

mg <- AddModuleScore_UCell(mg, features = gene_sets, ncores = 1)

# ---------------------------------------------------------------------------
# 2. Monocle3 拟时序
# ---------------------------------------------------------------------------
# 手工构建 cell_data_set（避免依赖 SeuratWrappers）
mg <- JoinLayers(mg)
counts <- GetAssayData(mg, assay = "RNA", layer = "counts")
cds <- new_cell_data_set(
  counts,
  cell_metadata = mg@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(mg),
                             row.names = rownames(mg))
)
# 把 Seurat 的 UMAP 坐标传给 monocle3
reducedDims(cds)[["UMAP"]] <- Embeddings(mg, "umap")

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

# 起点：LDAM_ref_UCell 最低分 5% 的细胞
scores <- mg@meta.data
root_cells <- rownames(scores)[scores$LDAM_ref_UCell <= quantile(scores$LDAM_ref_UCell, 0.05)]
cds <- order_cells(cds, root_cells = root_cells)

pt <- pseudotime(cds)
scores$pseudotime <- pt[rownames(scores)]
scores <- scores[!is.na(scores$pseudotime), ]
message("进入轨迹的细胞数: ", nrow(scores))

# ---------------------------------------------------------------------------
# 3. 轨迹动态：各评分沿拟时序的 Spearman 相关（单调性检验）
# ---------------------------------------------------------------------------
score_cols <- paste0(names(gene_sets), "_UCell")
traj_res <- bind_rows(lapply(score_cols, function(cc) {
  ct <- cor.test(scores$pseudotime, scores[[cc]], method = "spearman")
  tibble(module = str_remove(cc, "_UCell"), rho = unname(ct$estimate), P = ct$p.value)
}))
traj_res$FDR <- p.adjust(traj_res$P, "BH")
print(traj_res, n = Inf)
write_csv(traj_res, "results/31_pseudotime/31_trajectory_correlations.csv")

# 关键检验：LDAM 沿轨迹上升（定义性，应显著）；
#           铁死亡、Exosome_uptake、Cargo_loading 是否随之上升
# ---------------------------------------------------------------------------
# 4. 图
# ---------------------------------------------------------------------------
p_traj <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                     label_leaves = FALSE, label_branch_points = FALSE,
                     trajectory_graph_color = "grey30") +
  ggtitle("Microglia pseudotime (rooted at LDAM-low, GSE190452)")

dyn <- scores %>%
  select(pseudotime, all_of(score_cols)) %>%
  pivot_longer(-pseudotime, names_to = "module", values_to = "score") %>%
  mutate(module = str_remove(module, "_UCell"))

p_dyn <- ggplot(dyn, aes(pseudotime, score, color = module)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.6, linewidth = 0.9) +
  scale_color_manual(values = c("LDAM_ref" = "#D9A441", "Ferroptosis_ref" = "#C9828B",
                                "Exosome_uptake" = "#0E3A6B", "Cargo_loading" = "#4C956C",
                                "Rab_secretory" = "#8E7CC3")) +
  labs(title = "Module scores along the LDAM-low -> LDAM-high trajectory",
       x = "pseudotime (LDAM-low -> LDAM-high)", y = "UCell score")

pdf("results/31_pseudotime/31_pseudotime.pdf", width = 11, height = 5)
print(p_traj); print(p_dyn)
dev.off()

message("完成。判读：")
message("- LDAM_ref rho>0 为定义性阳性对照")
message("- Exosome_uptake / Cargo_loading rho>0 且 FDR<0.05 => 外泌体机器沿病程轴共同激活")
message("- 若铁死亡 rho 也>0，说明三者在同一演化轴上耦联推进")
scores <- mg@meta.data
mods <- intersect(paste0(c("Exosome_uptake","Cargo_loading","Rab_secretory",
                           "Syntenin_ectosome","ESCRT_core","SNARE_fusion",
                           "Ceramide_nSMase","Ferroptosis_ref"), "_UCell"),
                  colnames(scores))
cont_res <- bind_rows(lapply(mods, function(cc) {
  ct <- suppressWarnings(cor.test(scores$LDAM_ref_UCell, scores[[cc]], method = "spearman"))
  tibble(module = str_remove(cc, "_UCell"), rho = unname(ct$estimate), P = ct$p.value)
}))
cont_res$FDR <- p.adjust(cont_res$P, "BH")
print(cont_res, n = Inf)
write_csv(cont_res, "results/31_pseudotime/31_LDAM_continuous_correlations.csv")
