# ============================================================================
# Script      : 06_scRNA_GSE190452_main.R
# Purpose     : Core snRNA-seq analysis of GSE190452 (TLE-HS cortex, 4 vs 4): integration, annotation, LDAM/ferroptosis scoring, CellChat.
# Manuscript  : Methods 2.4 / Results 3.4-3.6
# Input       : GSE190452_RAW.tar (10x matrices); genesets/
# Output      : results/06_* (annotated Seurat rds, scores, CellChat); figures
# Dependencies: AUCell, CellChat, Seurat, SeuratWrappers, UCell, data.table, dplyr, ggplot2, ggpubr, harmony, monocle3, patchwork, stringr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================

# ==============================================================================
# 06_scRNA_GSE190452_main.R  【课题核心脚本】
# GSE190452：人 TLE（伴海马硬化）颞叶皮层 snRNA-seq，4 例 TLE vs 4 例非癫痫对照
# 修正版 v2：修复 cell_type 后缀、LDAM_level 转移、Seurat v5 JoinLayers
# ==============================================================================
options(repos = c(CRAN = "https://mirrors.westlake.edu.cn/cran"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
## Bioconductor 依赖
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# BiocManager::install(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats",
#                        "limma", "batchelor", "HDF5Array", "terra", "ggrastr"),
#                      ask = FALSE, update = FALSE)

## remotes + monocle3 + SeuratWrappers
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# install.packages("remotes")
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# remotes::install_github("cole-trapnell-lab/monocle3")
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# remotes::install_github("satijalab/seurat-wrappers")

library(monocle3); library(SeuratWrappers)
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results"); FIG_DIR <- file.path(BASE_DIR, "figures")
GS_DIR   <- file.path(BASE_DIR, "genesets")

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(dplyr); library(ggplot2)
  library(data.table); library(stringr); library(UCell); library(AUCell)
  library(patchwork); library(ggpubr)
})
set.seed(2026)

LDAM_genes <- readLines(file.path(GS_DIR, "LDAM_geneset.txt"))
Ferroptosis_genes <- readLines(file.path(GS_DIR, "Ferroptosis_geneset.txt"))

## ---------- 6.1 解压并读取 8 个样本的 10x 矩阵 ----------
untar_dir <- file.path(BASE_DIR, "GSE190452")
if (!dir.exists(untar_dir)) {
  dir.create(untar_dir)
  untar(file.path(BASE_DIR, "GSE190452_RAW.tar"), exdir = untar_dir)
}
all_files <- list.files(untar_dir, recursive = TRUE)
print(all_files)

matrix_files <- list.files(untar_dir, pattern = "matrix\\.mtx(\\.gz)?$",
                           recursive = TRUE, full.names = TRUE)
matrix_files <- matrix_files[!grepl("std_", basename(dirname(matrix_files)))]
matrix_files <- list.files(untar_dir, pattern = "matrix\\.mtx(\\.gz)?$",
                           recursive = TRUE, full.names = TRUE)
matrix_files <- matrix_files[!grepl("std_", basename(dirname(matrix_files)))]
print(matrix_files)        # 应该正好 8 个，且不包含 std_

seurat_list <- lapply(matrix_files, function(mf) {
  m <- read_10x_flexible(mf)
  sname <- sub(".*(GSM\\d+[^/\\\\]*).*", "\\1", mf)
  CreateSeuratObject(m, project = sname, min.cells = 3, min.features = 200)
})
names(seurat_list) <- sapply(seurat_list, function(x) x$orig.ident[1])
# [debug] print(names(seurat_list))
read_10x_flexible <- function(mtx_file) {
  d <- dirname(mtx_file)
  tmp <- file.path(d, paste0("std_", basename(d), "_", tools::file_path_sans_ext(basename(mtx_file))))
  dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
  prefix <- sub("_?matrix\\.mtx(\\.gz)?$", "", basename(mtx_file))
  for (suffix in c("barcodes.tsv", "features.tsv", "genes.tsv", "matrix.mtx")) {
    cand <- list.files(d, pattern = paste0(suffix, "(\\.gz)?$"), full.names = TRUE)
    cand <- cand[grepl(prefix, basename(cand), fixed = TRUE)]
    if (length(cand) == 0 && suffix == "features.tsv")
      cand <- list.files(d, pattern = "genes\\.tsv(\\.gz)?$", full.names = TRUE)
    if (length(cand) > 0)
      file.copy(cand[1], file.path(tmp, paste0(suffix, ifelse(grepl("\\.gz$", cand[1]), ".gz", ""))),
                overwrite = TRUE)
  }
  Read10X(tmp)
}

seurat_list <- lapply(matrix_files, function(mf) {
  m <- read_10x_flexible(mf)
  sname <- sub(".*(GSM\\d+[^/\\\\]*).*", "\\1", mf)
  CreateSeuratObject(m, project = sname, min.cells = 3, min.features = 200)
})
names(seurat_list) <- sapply(seurat_list, function(x) x$orig.ident[1])
# [debug] print(names(seurat_list))

for (i in seq_along(seurat_list)) {
  nm <- names(seurat_list)[i]
  seurat_list[[i]]$group <- ifelse(str_detect(tolower(nm), "normal|control"), "Control", "TLE")
}
print(sapply(seurat_list, function(x) table(x$group)))

## ---------- 6.2 合并 + 质控 ----------
obj <- merge(seurat_list[[1]], seurat_list[-1])
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 10)

p_qc <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "orig.ident", ncol = 3, pt.size = 0)
ggsave(file.path(FIG_DIR, "06_QC_violin.pdf"), p_qc, width = 12, height = 5)

## ---------- 6.3 标准化 + Harmony 去批次 + 聚类 ----------
obj <- NormalizeData(obj) |>
  FindVariableFeatures(nfeatures = 3000) |>
  ScaleData(vars.to.regress = c("percent.mt")) |>
  RunPCA(npcs = 50)
obj <- RunHarmony(obj, group.by.vars = "orig.ident")
obj <- obj |>
  RunUMAP(reduction = "harmony", dims = 1:30) |>
  FindNeighbors(reduction = "harmony", dims = 1:30) |>
  FindClusters(resolution = 0.5)

p_umap <- DimPlot(obj, group.by = "orig.ident") + DimPlot(obj, label = TRUE)
ggsave(file.path(FIG_DIR, "06_UMAP_overview.pdf"), p_umap, width = 14, height = 6)

## ---------- 6.4 细胞类型注释 ----------
markers <- list(
  Excitatory_neuron = c("SLC17A7", "CAMK2A", "SATB2"),
  Inhibitory_neuron = c("GAD1", "GAD2", "SLC32A1"),
  Astrocyte         = c("AQP4", "GFAP", "SLC1A2"),
  Microglia         = c("P2RY12", "TMEM119", "CX3CR1", "CSF1R"),
  Oligodendrocyte   = c("MBP", "PLP1", "MOG"),
  OPC               = c("VCAN", "PDGFRA", "CSPG4"),
  Endothelial       = c("CLDN5", "VWF", "FLT1")
)
for (ct in names(markers)) {
  g <- intersect(markers[[ct]], rownames(obj))
  if (length(g) >= 2) obj <- AddModuleScore(obj, features = list(g), name = ct)
}
score_cols <- paste0(names(markers), "1")
score_cols <- score_cols[score_cols %in% colnames(obj@meta.data)]
avg_scores <- obj@meta.data |>
  group_by(seurat_clusters) |>
  summarise(across(all_of(score_cols), mean))
print(avg_scores)

cluster_anno <- avg_scores |>
  tidyr::pivot_longer(-seurat_clusters, names_to = "ct", values_to = "score") |>
  group_by(seurat_clusters) |> slice_max(score, n = 1) |>
  select(seurat_clusters, cell_type = ct)
obj$cell_type <- cluster_anno$cell_type[match(obj$seurat_clusters, cluster_anno$seurat_clusters)]
obj$cell_type <- sub("\\d+$", "", obj$cell_type)        # 【修复1】去掉 "1" 后缀
print(table(obj$cell_type))

p_anno <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("GSE190452 cell type annotation")
ggsave(file.path(FIG_DIR, "06_UMAP_celltype.pdf"), p_anno, width = 8, height = 6)

## ---------- 6.5 【Go/No-Go 检查点 1】小胶质细胞提取 ----------
mg <- subset(obj, subset = cell_type == "Microglia")
n_mg <- ncol(mg)
print(table(mg$group))
message("【Go/No-Go 检查点 1】小胶质细胞核总数 = ", n_mg,
        "（方案阈值 >= 1000）", ifelse(n_mg >= 1000, "→ 达标 GO", "→ 未达标，需评估"))

## ---------- 6.6 LDAM / 铁死亡评分（UCell） ----------
mg <- AddModuleScore_UCell(mg, features = list(LDAM = LDAM_genes,
                                               Ferroptosis = Ferroptosis_genes))
# 【Go/No-Go 检查点 2】细胞水平 Wilcoxon + 样本水平 pseudobulk 双验证
wil <- wilcox.test(LDAM_UCell ~ group, data = mg@meta.data)
pb <- mg@meta.data |>
  group_by(orig.ident, group) |>
  summarise(LDAM = mean(LDAM_UCell), .groups = "drop")
tt <- t.test(LDAM ~ group, data = pb)
message("【Go/No-Go 检查点 2】LDAM 评分 TLE vs 对照：细胞水平 Wilcoxon P = ",
        signif(wil$p.value, 3), "；样本水平 t 检验 P = ", signif(tt$p.value, 3))
print(tapply(pb$LDAM, pb$group, mean))     # 方向：TLE 应 > Control
fwrite(pb, file.path(OUT_DIR, "06_microglia_LDam_score_per_sample.csv"))

p_score <- ggplot(mg@meta.data, aes(group, LDAM_UCell, fill = group)) +
  geom_violin(scale = "width") + geom_boxplot(width = 0.15, outlier.size = 0.3) +
  scale_fill_manual(values = c("Control" = "#7BA7BC", "TLE" = "#C0392B")) +
  labs(title = paste0("Microglia LDAM score (UCell)\nsample-level P = ",
                      signif(tt$p.value, 3))) + theme_bw()
ggsave(file.path(FIG_DIR, "06_microglia_LDam_score.pdf"), p_score, width = 6, height = 6)

p_pb <- ggplot(pb, aes(group, LDAM, fill = group)) +
  geom_boxplot(width = 0.4) + geom_point(size = 3) +
  scale_fill_manual(values = c("Control" = "#7BA7BC", "TLE" = "#C0392B")) +
  labs(title = "Sample-level LDAM score (pseudobulk)") + theme_bw()
ggsave(file.path(FIG_DIR, "06_microglia_LDam_score_pseudobulk.pdf"), p_pb, width = 5, height = 6)

## ---------- 6.7 小胶质细胞重聚类 + LDAM-high/low 亚群 ----------
mg <- NormalizeData(mg) |> FindVariableFeatures(nfeatures = 2000) |>
  ScaleData() |> RunPCA(npcs = 30)
mg <- RunHarmony(mg, group.by.vars = "orig.ident")
mg <- mg |> RunUMAP(reduction = "harmony", dims = 1:20) |>
  FindNeighbors(reduction = "harmony", dims = 1:20) |>
  FindClusters(resolution = 0.3)

mg$LDAM_level <- ifelse(mg$LDAM_UCell > median(mg$LDAM_UCell), "LDAM-high", "LDAM-low")

# 亚群占比：LDAM-high 比例 TLE vs 对照（样本水平统计）
prop_tab <- prop.table(table(mg$orig.ident, mg$LDAM_level), margin = 1) |>
  as.data.frame()
colnames(prop_tab) <- c("sample", "LDAM_level", "proportion")
prop_tab$group <- unique(mg@meta.data[, c("orig.ident", "group")])$group[
  match(prop_tab$sample, unique(mg@meta.data[, c("orig.ident", "group")])$orig.ident)]
fwrite(prop_tab, file.path(OUT_DIR, "06_LDAm_high_low_proportion_per_sample.csv"))
prop_hi <- prop_tab[prop_tab$LDAM_level == "LDAM-high", ]
wil_prop <- wilcox.test(proportion ~ group, data = prop_hi)
message("LDAM-high 小胶质细胞占比 TLE vs 对照：Wilcoxon P = ", signif(wil_prop$p.value, 3))
print(prop_hi)

p_prop <- ggplot(prop_hi, aes(group, proportion, fill = group)) +
  geom_boxplot(width = 0.4) + geom_point(size = 3) +
  scale_fill_manual(values = c("Control" = "#7BA7BC", "TLE" = "#C0392B")) +
  labs(title = "LDAM-high microglia proportion") + theme_bw()
ggsave(file.path(FIG_DIR, "06_LDam_high_low_proportion.pdf"), p_prop, width = 5, height = 6)

p_mg_umap <- DimPlot(mg, group.by = "LDAM_level",
                     cols = c("LDAM-high" = "#C0392B", "LDAM-low" = "#7BA7BC")) +
  FeaturePlot(mg, features = "LDAM_UCell") & theme_bw()
ggsave(file.path(FIG_DIR, "06_microglia_UMAP_LDam.pdf"), p_mg_umap, width = 14, height = 6)

## ---------- 6.8 LDAM vs 铁死亡相关性（单细胞层面） ----------
cor_cell <- cor.test(mg$LDAM_UCell, mg$Ferroptosis_UCell, method = "spearman")
message("单细胞层面 LDAM-铁死亡相关：rho = ", round(cor_cell$estimate, 3),
        "，P = ", signif(cor_cell$p.value, 3))

p_corr <- ggplot(mg@meta.data, aes(LDAM_UCell, Ferroptosis_UCell, color = group)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("Control" = "#7BA7BC", "TLE" = "#C0392B")) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("Spearman r = ", round(cor_cell$estimate, 3),
                          ", P = ", signif(cor_cell$p.value, 3))) +
  labs(title = "Microglia LDAM vs Ferroptosis score") + theme_bw()
ggsave(file.path(FIG_DIR, "06_LDam_Ferroptosis_correlation_sc.pdf"), p_corr, width = 7, height = 6)

p_ferr <- ggplot(mg@meta.data, aes(LDAM_level, Ferroptosis_UCell, fill = LDAM_level)) +
  geom_violin(scale = "width") + geom_boxplot(width = 0.15, outlier.size = 0.3) +
  stat_compare_means(method = "wilcox.test") +
  scale_fill_manual(values = c("LDAM-high" = "#C0392B", "LDAM-low" = "#7BA7BC")) +
  labs(title = "Ferroptosis score: LDAM-high vs LDAM-low") + theme_bw()
ggsave(file.path(FIG_DIR, "06_ferroptosis_LDamHighLow.pdf"), p_ferr, width = 6, height = 6)

## ---------- 6.9 hub 基因在小胶质细胞中的表达 ----------
hub_file <- file.path(OUT_DIR, "04_hub_genes.txt")
if (file.exists(hub_file)) {
  hub_genes <- intersect(readLines(hub_file), rownames(mg))
  if (length(hub_genes) > 0) {
    p_hub <- VlnPlot(mg, features = hub_genes, group.by = "group",
                     pt.size = 0, ncol = min(3, length(hub_genes)))
    ggsave(file.path(FIG_DIR, "06_hub_genes_violin.pdf"), p_hub,
           width = 4 * min(3, length(hub_genes)), height = 4)
  }
}

saveRDS(mg,  file.path(OUT_DIR, "06_GSE190452_microglia.rds"))
saveRDS(obj, file.path(OUT_DIR, "06_GSE190452_annotated_full.rds"))

## ---------- 6.10 Monocle3 拟时序 ----------
if (requireNamespace("monocle3", quietly = TRUE) &
    requireNamespace("SeuratWrappers", quietly = TRUE)) {
  library(monocle3); library(SeuratWrappers)
  cds <- as.cell_data_set(mg)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  root_cells <- colnames(mg)[order(mg$LDAM_UCell)[1:50]]
  cds <- order_cells(cds, root_cells = root_cells)
  p_pseudo <- plot_cells(cds, color_cells_by = "pseudotime",
                         label_cell_groups = FALSE, label_leaves = FALSE,
                         label_branch_points = FALSE) +
    ggtitle("Microglia pseudotime (homeostatic -> LDAM)")
  ggsave(file.path(FIG_DIR, "06_monocle3_pseudotime.pdf"), p_pseudo, width = 7, height = 6)
  ptime <- pseudotime(cds)
  df_pt <- data.frame(pseudotime = ptime, LDAM = mg$LDAM_UCell,
                      Ferroptosis = mg$Ferroptosis_UCell)
  p_pt_line <- ggplot(df_pt, aes(pseudotime, LDAM)) +
    geom_point(alpha = 0.3, size = 0.5) + geom_smooth(color = "#C0392B") +
    geom_smooth(aes(y = Ferroptosis), color = "#2C5F8A") +
    labs(title = "LDAM / Ferroptosis score along pseudotime") + theme_bw()
  ggsave(file.path(FIG_DIR, "06_pseudotime_score_trend.pdf"), p_pt_line, width = 7, height = 5)
} else {
  message("monocle3 未安装，跳过拟时序分析")
}

## ---------- 6.11 CellChat 细胞通讯 ----------
if (requireNamespace("CellChat", quietly = TRUE)) {
  library(CellChat)
  obj <- JoinLayers(obj)                                   # 【修复3】合并 v5 多层
  obj$LDAM_level <- NA_character_                          # 【修复2】转回完整对象
  obj$LDAM_level[colnames(mg)] <- mg$LDAM_level
  obj$chat_group <- ifelse(obj$cell_type == "Microglia",
                           ifelse(obj$LDAM_level %in% "LDAM-high",
                                  "Microglia_LDamHigh", "Microglia_LDamLow"),
                           obj$cell_type)
  obj$chat_group[obj$cell_type == "Microglia" & is.na(obj$chat_group)] <- "Microglia_LDamLow"
  print(table(obj$chat_group))

  chat <- createCellChat(object = GetAssayData(obj, layer = "data"),
                         meta = obj@meta.data, group.by = "chat_group")
  chat@DB <- CellChatDB.human
  chat <- subsetData(chat) |> identifyOverExpressedGenes() |>
    identifyOverExpressedInteractions() |>
    computeCommunProb() |> filterCommunication(min.cells = 10) |>
    computeCommunProbPathway() |> aggregateNet()

  pdf(file.path(FIG_DIR, "06_CellChat_circle.pdf"), width = 8, height = 8)
  print(netVisual_circle(chat@net$count, vertex.weight = as.numeric(table(chat@idents)),
                         weight.scale = TRUE, label.edge = FALSE,
                         title.name = "Cell-cell communication (counts)"))
  dev.off()

  message("显著信号通路："); print(chat@netP$pathways)
  for (pw in c("IL1", "TNF", "COMPLEMENT")) {
    if (pw %in% chat@netP$pathways) {
      pdf(file.path(FIG_DIR, paste0("06_CellChat_", pw, ".pdf")), width = 8, height = 8)
      print(netVisual_aggregate(chat, signaling = pw, layout = "circle"))
      dev.off()
    }
  }
  saveRDS(chat, file.path(OUT_DIR, "06_CellChat_object.rds"))
} else {
  message("CellChat 未安装，跳过细胞通讯分析")
}

message("=== 06_scRNA_GSE190452_main.R 运行完毕 ===")
library(CellChat)
chat <- readRDS(file.path(OUT_DIR, "06_CellChat_object.rds"))
for (pw in c("SPP1", "Cholesterol", "TGFb")) {
  if (pw %in% chat@netP$pathways) {
    pdf(file.path(FIG_DIR, paste0("06_CellChat_", pw, ".pdf")), width = 8, height = 8)
    print(netVisual_aggregate(chat, signaling = pw, layout = "circle"))
    dev.off()
    cat(pw, "图已保存\n")
  }
}
