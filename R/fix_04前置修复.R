# ============================================================================
# Script      : fix_04前置修复.R
# Purpose     : LEGACY session-continuation helper (Harmony -> clustering -> UMAP -> annotation) from an earlier project; requires objects already in session.
# Manuscript  : not used in manuscript (legacy helper)
# Input       : existing session objects (merged Seurat after PCA)
# Output      : rds/figures/tables under a legacy WORK_DIR
# Dependencies: Matrix, Seurat, SeuratObject, ggplot2
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# Status      : LEGACY / DEPRECATED - kept for transparency, not part of the manuscript pipeline
# ============================================================================

## ============================================================
## 03 后半段完整版：Harmony → 聚类 → UMAP → 注释 → 保存
## 前提：会话中已有 merged，且已完成 PCA
## ============================================================

## ---- 0. 路径与依赖 ----
## [EN] Step 0: Paths and dependencies
# [legacy] original path pointed to an earlier glioma project; parameterised below
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root

if (!exists("DATA_DIR")) DATA_DIR <- Sys.getenv("GLIOMA_GEO_DIR", file.path(BASE_DIR, "data/external/glioma/GEO"))
if (!exists("WORK_DIR")) WORK_DIR <- file.path(DATA_DIR, "analysis")
RDS_DIR <- file.path(WORK_DIR, "rds")
FIG_DIR <- file.path(WORK_DIR, "figures")
TBL_DIR <- file.path(WORK_DIR, "tables")
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("FNN", quietly = TRUE))
  # [setup] one-time package installation commented out for submission;
  # [setup] run 00_setup_genesets.R (or install manually) before first use:
#   install.packages("FNN", repos = "https://mirrors.westlake.edu.cn/CRAN")

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(ggplot2); library(Matrix)
})
if (!exists("NPCS")) NPCS <- 30
stopifnot(exists("merged"))
message("细胞数: ", ncol(merged), " | reductions: ",
        paste(names(merged@reductions), collapse = ", "))

## ---- 1. Harmony 批次校正（已跑过则跳过）----
if ("harmony" %in% names(merged@reductions)) {
  message("已有 harmony reduction，跳过。")
} else {
  batch_var <- if ("dataset" %in% colnames(merged@meta.data)) "dataset" else "sample_id"
  message("批次变量: ", batch_var)

  pca_emb <- Embeddings(merged, "pca")[, 1:NPCS]
  meta_use <- merged@meta.data
  ok_cells <- rownames(pca_emb)[stats::complete.cases(pca_emb) &
                                  !is.na(meta_use[rownames(pca_emb), batch_var])]
  message("有效细胞: ", length(ok_cells), " / ", ncol(merged))
  pca_emb <- pca_emb[ok_cells, ]
  meta_use <- meta_use[ok_cells, , drop = FALSE]
  meta_use[[batch_var]] <- as.factor(meta_use[[batch_var]])

  hm_mat <- NULL
  if (requireNamespace("harmony", quietly = TRUE)) {
    message("运行 harmony（矩阵接口，自动适配版本）...")
    rh_args <- formals(harmony::RunHarmony)
    if ("group.by.vars" %in% names(rh_args)) {          # harmony >= 2.0
      hm_obj <- harmony::RunHarmony(pca_emb, meta_use,
                                    group.by.vars = batch_var, verbose = TRUE)
    } else {                                            # harmony 1.x
      hm_obj <- harmony::RunHarmony(pca_emb, meta_use,
                                    vars_use = batch_var, verbose = TRUE)
    }
    hm_mat <- if (is.matrix(hm_obj)) hm_obj else hm_obj$harmony_embeddings
    rm(hm_obj); gc(verbose = FALSE)
  }

  if (is.null(hm_mat)) {
    message("harmony 失败，回退到 Seurat::IntegrateLayers ...")
    merged <- IntegrateLayers(merged, method = HarmonyIntegration,
                              orig.reduction = "pca",
                              new.reduction = "harmony", verbose = FALSE)
  } else {
    ## 按位置回填（完整细胞数，NA 行已被剔除的需对齐）
    full_hm <- matrix(NA_real_, nrow = ncol(merged), ncol = ncol(hm_mat),
                      dimnames = list(colnames(merged), colnames(hm_mat)))
    idx_ok <- match(rownames(hm_mat), colnames(merged))
    full_hm[idx_ok, ] <- hm_mat
    merged[["harmony"]] <- CreateDimReducObject(
      embeddings = full_hm, key = "harmony_", assay = "RNA")
    rm(hm_mat, full_hm, pca_emb, meta_use); gc(verbose = FALSE)
  }
  message("harmony 完成。dims: ", ncol(Embeddings(merged, "harmony")))
}

## ---- 2. 子采样聚类 + UMAP ----
message("\n--- 子采样聚类与 UMAP ---")
set.seed(123)
N_SKETCH <- 50000
idx_sketch <- if (ncol(merged) > N_SKETCH) sample(ncol(merged), N_SKETCH) else seq_len(ncol(merged))
sketch_cells <- colnames(merged)[idx_sketch]

sub <- DietSeurat(merged[, sketch_cells],
                  assays = "RNA", dimreducs = "harmony", graphs = character())
sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:NPCS, verbose = FALSE)
sub <- FindClusters(sub, resolution = 0.8, verbose = FALSE)
sub <- RunUMAP(sub, reduction = "harmony", dims = 1:NPCS,
               n.neighbors = 30, min.dist = 0.3, verbose = FALSE)
message("子集聚类: ", length(levels(sub$seurat_clusters)), " 个簇")

## ---- 3. 全量 UMAP 投影 + 簇标签回填（全位置索引）----
idx_sub <- match(colnames(sub), colnames(merged))
if (any(is.na(idx_sub))) stop("子集细胞位置匹配失败，请把输出贴给我")

full_umap <- matrix(0, nrow = ncol(merged), ncol = 2,
                    dimnames = list(colnames(merged), c("umap_1", "umap_2")))
full_umap[idx_sub, ] <- as.matrix(Embeddings(sub, "umap"))

merged$seurat_clusters <- NA_character_
merged$seurat_clusters[idx_sub] <- as.character(sub$seurat_clusters)

rest <- setdiff(seq_len(ncol(merged)), idx_sub)
chunk <- 20000
if (length(rest) > 0) {
  cl_ref   <- as.character(sub$seurat_clusters)
  sub_harm <- Embeddings(sub, "harmony")[, 1:NPCS]
  sub_umap <- Embeddings(sub, "umap")
  message("投影其余 ", length(rest), " 个细胞...")
  for (i in seq(1, length(rest), by = chunk)) {
    cc <- rest[i:min(i + chunk - 1, length(rest))]
    kn <- FNN::get.knnx(sub_harm, Embeddings(merged, "harmony")[cc, 1:NPCS], k = 1)
    nn <- kn$nn.index[, 1]
    full_umap[cc, ] <- sub_umap[nn, ]
    merged$seurat_clusters[cc] <- cl_ref[nn]
    if (i %% (chunk * 5) == 1) message("  已处理 ", min(i + chunk - 1, length(rest)),
                                       " / ", length(rest))
    gc(verbose = FALSE)
  }
}
merged$seurat_clusters <- factor(merged$seurat_clusters)
message("全量簇标签完成: ", length(levels(merged$seurat_clusters)), " 个簇")

merged[["umap"]] <- CreateDimReducObject(
  embeddings = full_umap, key = "umap_", assay = "RNA")
rm(sub, full_umap); gc(verbose = FALSE)

## ---- 4. 模块打分注释（逐层稀疏，不 JoinLayers）----
message("\n--- 模块打分注释 ---")
marker_sets <- list(
  Myeloid     = c("PTPRC","CD68","AIF1","ITGAM","FCGR3A","CD163","MRC1","LYZ","CSF1R"),
  Microglia   = c("P2RY12","TMEM119","SALL1","CX3CR1","GPR34","HEXB","SELPLG"),
  Mono_Macro  = c("CD14","FCN1","VCAN","S100A8","S100A9","CCR2","LILRB2"),
  T_cell      = c("CD3D","CD3E","CD8A","CD4","IL7R","TRAC"),
  B_cell      = c("CD79A","MS4A1","CD19","MZB1","IGKC"),
  NK          = c("NKG7","GNLY","KLRD1","NCAM1"),
  Tumor       = c("EGFR","PDGFRA","SOX2","NES","VIM","CDK4","OLIG1","OLIG2","ASCL1"),
  Oligodendro = c("MBP","MOG","PLP1","MAG","CLDN11"),
  Astrocyte   = c("GFAP","AQP4","ALDH1L1","SLC1A3"),
  Endothelial = c("PECAM1","VWF","KDR","CLDN5","FLT1"),
  Pericyte    = c("RGS5","PDGFRB","ACTA2","MCAM"),
  Fibroblast  = c("COL1A1","DCN","COL3A1","LUM")
)

score_module <- function(obj, genes) {
  layers <- SeuratObject::Layers(obj[["RNA"]], search = "^data")
  if (length(layers) == 0) stop("没有 data 层，请先对每个数据集 NormalizeData")
  acc <- numeric(ncol(obj)); names(acc) <- colnames(obj)
  for (ly in layers) {
    mat <- LayerData(obj, assay = "RNA", layer = ly)
    g_use <- intersect(rownames(mat), genes)
    if (length(g_use) < 2) next
    acc[colnames(mat)] <- Matrix::colMeans(mat[g_use, , drop = FALSE])
  }
  acc
}

for (nm in names(marker_sets)) {
  merged[[paste0("ms_", nm)]] <- score_module(merged, marker_sets[[nm]])
  message("  打分完成: ", nm)
}

ms_mat <- as.matrix(merged@meta.data[, paste0("ms_", names(marker_sets))])
ms_mat[, "ms_Myeloid"] <- pmax(ms_mat[, "ms_Microglia"], ms_mat[, "ms_Mono_Macro"])
merged$marker_based_anno <- factor(names(marker_sets)[max.col(ms_mat)],
                                   levels = names(marker_sets))
message("\n注释分布:")
print(table(merged$marker_based_anno))

## ---- 5. 出图 ----
message("\n--- 出图 ---")
DimPlot(merged, reduction = "umap", group.by = "dataset",
        raster = TRUE, shuffle = TRUE) + ggtitle("Batch check (by dataset)")
ggsave(file.path(FIG_DIR, "03_umap_by_dataset.png"), width = 9, height = 7, dpi = 300)

DimPlot(merged, reduction = "umap", group.by = "seurat_clusters",
        raster = TRUE, label = TRUE) + ggtitle("Clusters")
ggsave(file.path(FIG_DIR, "03_umap_clusters.png"), width = 9, height = 7, dpi = 300)

DimPlot(merged, reduction = "umap", group.by = "marker_based_anno",
        raster = TRUE, label = TRUE) + ggtitle("Marker-based annotation")
ggsave(file.path(FIG_DIR, "03_umap_annotation.png"), width = 10, height = 7, dpi = 300)

avg <- aggregate(merged@meta.data[, c("ms_Myeloid","ms_Microglia","ms_Mono_Macro",
                                      "ms_Tumor","ms_T_cell")],
                 by = list(cluster = merged$seurat_clusters), FUN = mean)
data.table::fwrite(avg, file.path(TBL_DIR, "03_cluster_module_scores.csv"))

## ---- 6. 保存 ----
saveRDS(merged, file.path(RDS_DIR, "merged_annotated.rds"))
message("\n=== 03 完成 ===")
message("保存: ", file.path(RDS_DIR, "merged_annotated.rds"))
message("reductions: ", paste(names(merged@reductions), collapse = ", "))
message("簇数: ", length(levels(merged$seurat_clusters)))
message("髓系细胞数(粗注释): ",
        sum(merged$marker_based_anno %in% c("Myeloid","Microglia","Mono_Macro")))
