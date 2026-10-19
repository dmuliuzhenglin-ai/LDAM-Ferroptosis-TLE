# ============================================================================
# Script      : 00_setup_genesets.R
# Purpose     : One-time setup: check/install dependencies and build + freeze the LDAM and ferroptosis gene sets.
# Manuscript  : Methods 2.1 (Datasets and gene sets)
# Input       : MSigDB collections via msigdbr; built-in LDAM gene list
# Output      : genesets/ (frozen gene-set rds/csv for archival)
# Dependencies: data.table, dplyr, ggplot2, msigdbr, stringr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

# ==============================================================================
# 00_setup_genesets.R
# 课题：LDAM 联合铁死亡在海马起源颞叶癫痫（HTLE）中的生信研究（方案 v1.3）
# 功能：① 检查/安装全部依赖包 ② 构建并固定保存 LDAM 与铁死亡基因集
# 运行顺序：第 1 个运行，只需运行一次
# ==============================================================================

## ---------- 0.1 全局路径设置（所有脚本共用同一套路径） ----------
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
DATA_DIR <- BASE_DIR                                   # 数据文件直接在项目文件夹中
OUT_DIR  <- file.path(BASE_DIR, "results")             # 结果表格
FIG_DIR  <- file.path(BASE_DIR, "figures")             # 图片
GS_DIR   <- file.path(BASE_DIR, "genesets")            # 固定基因集（论文存档用）
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GS_DIR,  showWarnings = FALSE, recursive = TRUE)

## ---------- 0.2 依赖包检查与安装 ----------
cran_pkgs <- c("data.table", "dplyr", "tidyr", "stringr", "ggplot2", "ggrepel",
               "pheatmap", "patchwork", "ggpubr", "randomForest", "pROC",
               "glmnet", "e1071", "caret", "msigdbr", "WGCNA", "Seurat",
               "harmony", "Matrix", "cowplot", "RColorBrewer")
bioc_pkgs <- c("limma", "DESeq2", "sva", "clusterProfiler", "GSVA", "GSEABase",
               "org.Hs.eg.db", "GEOquery", "AUCell", "UCell", "biomaRt",
               "enrichplot", "impute")
other_pkgs <- c("monocle3", "SeuratWrappers", "CellChat")   # 见下方特殊安装说明

install_if_missing <- function(pkgs, installer) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("安装缺失包: ", p)
      tryCatch(installer(p), error = function(e)
        message("!! ", p, " 自动安装失败，请手动安装: ", conditionMessage(e)))
    }
  }
}
install_if_missing(cran_pkgs, install.packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install_if_missing(bioc_pkgs, BiocManager::install)

# monocle3 / SeuratWrappers / CellChat 为 GitHub 或 Bioc 特殊源：
# if (!requireNamespace("monocle3", quietly = TRUE))
#   BiocManager::install(c('BiocGenerics','DelayedArray','DelayedMatrixStats',
#                          'limma','S4Vectors','SingleCellExperiment',
#                          'SummarizedExperiment','batchelor','HDF5Array',
#                          'terra','ggrastr'))
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# if (!requireNamespace("monocle3", quietly = TRUE)) devtools::install_github('cole-trapnell-lab/monocle3')
# if (!requireNamespace("SeuratWrappers", quietly = TRUE)) remotes::install_github('satijalab/seurat-wrappers')
# if (!requireNamespace("CellChat", quietly = TRUE)) devtools::install_github("jinworks/CellChat")

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(ggplot2)
})

## ---------- 0.3 构建 LDAM 基因集（MSigDB 官方名称 + 文献基因） ----------
# 方案 v1.3 已核实：MSigDB 中正确名称为 GOCC_LIPID_DROPLET /
# GOBP_TRIGLYCERIDE_METABOLIC_PROCESS / HALLMARK_FATTY_ACID_METABOLISM
library(msigdbr)

get_msig <- function(category, subcategory, gs) {
  msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory) |>
    filter(gs_name == gs) |> pull(gene_symbol) |> unique()
}

msig_lipid_droplet <- get_msig("C5", "GO:CC", "GOCC_LIPID_DROPLET")                       # 102 基因
msig_tg_process    <- get_msig("C5", "GO:BP", "GOBP_TRIGLYCERIDE_METABOLIC_PROCESS")      # 119 基因
msig_fa_hallmark   <- get_msig("H",  "",      "HALLMARK_FATTY_ACID_METABOLISM")           # 158 基因
message("GOCC_LIPID_DROPLET: ", length(msig_lipid_droplet),
        " | GOBP_TRIGLYCERIDE_METABOLIC_PROCESS: ", length(msig_tg_process),
        " | HALLMARK_FATTY_ACID_METABOLISM: ", length(msig_fa_hallmark))

# 文献报道 LDAM 标志基因（方案 2.4 节，存档备查）
ldam_literature <- c("APOE","TREM2","CD9","PLIN2","PLIN3","ACSL1","ACSL4","CPT1A",
                     "DGAT1","DGAT2","LPCAT2","LPL","CH25H","FABP5",
                     "CTSD","LIPA","SPP1","GPNMB","CSTD","TYROBP","C1QA","C1QB","C1QC")

LDAM_genes <- unique(c(msig_lipid_droplet, msig_tg_process, msig_fa_hallmark, ldam_literature))
message("LDAM 基因集（并集去重）: ", length(LDAM_genes), " 个基因")

## ---------- 0.4 构建铁死亡基因集（FerrDb V2/V3 + MSigDB） ----------
# FerrDb 需手动下载：打开 http://www.zhounan.org/ferrdb/current/operations/download.html
# 下载 driver / suppressor / marker 的 csv，放入 genesets/ferrdb/ 目录
ferr_dir <- file.path(GS_DIR, "ferrdb"); dir.create(ferr_dir, showWarnings = FALSE)
ferr_files <- list.files(ferr_dir, pattern = "\\.csv$", full.names = TRUE)

ferrdb_genes <- character(0)
if (length(ferr_files) > 0) {
  for (f in ferr_files) {
    tmp <- tryCatch(fread(f), error = function(e) NULL)
    if (!is.null(tmp)) {
      # FerrDb 表格的基因列通常名为 symbol / gene_symbol，做兼容处理
      gcol <- grep("symbol|gene", names(tmp), ignore.case = TRUE, value = TRUE)[1]
      if (!is.na(gcol)) ferrdb_genes <- c(ferrdb_genes, toupper(tmp[[gcol]]))
    }
  }
  ferrdb_genes <- unique(ferrdb_genes)
  message("FerrDb 本地文件读取: ", length(ferrdb_genes), " 个基因")
} else {
  message("未发现 FerrDb csv（genesets/ferrdb/ 为空），回退使用 MSigDB 铁死亡通路基因集")
}

msig_ferroptosis <- unique(c(
  get_msig("C5", "GO:BP", "GOBP_FERROPTOSIS"),
  tryCatch(get_msig("C2", "CP:WIKIPATHWAYS", "WP_FERROPTOSIS"), error = function(e) character(0)),
  get_msig("C2", "CP:REACTOME", "REACTOME_FERROPTOSIS")
))
# 铁死亡核心基因兜底（FerrDb 经典分子，保证网络不佳时脚本仍可运行）
ferroptosis_core <- c("GPX4","SLC7A11","SLC3A2","ACSL4","LPCAT3","ALOX15","ALOX12",
                      "FTH1","FTL","TFRC","NCOA4","HMOX1","NFE2L2","KEAP1","GCLC",
                      "GCLM","SAT1","CBS","NOX1","NOX4","PTGS2","CHAC1","ATF4","DDIT3",
                      "VDAC2","VDAC3","GOT1","GLS2","CARS1","FSP1","AIFM2","DHODH",
                      "GCH1","NFS1","ISCU","ABCC1","MT1G","HSPB1","BAP1","BRD4")

Ferroptosis_genes <- unique(c(ferrdb_genes, msig_ferroptosis, ferroptosis_core))
message("铁死亡基因集（并集去重）: ", length(Ferroptosis_genes), " 个基因")

## ---------- 0.5 固定保存基因集（全文统一调用，论文存档备查） ----------
writeLines(LDAM_genes,        file.path(GS_DIR, "LDAM_geneset.txt"))
writeLines(Ferroptosis_genes, file.path(GS_DIR, "Ferroptosis_geneset.txt"))
union_genes <- unique(c(LDAM_genes, Ferroptosis_genes))
writeLines(union_genes, file.path(GS_DIR, "LDAM_Ferroptosis_union.txt"))

fwrite(data.frame(gene = LDAM_genes,        set = "LDAM"),
       file.path(GS_DIR, "geneset_archive_long.csv"))
fwrite(data.frame(gene = Ferroptosis_genes, set = "Ferroptosis"),
       file.path(GS_DIR, "geneset_archive_long.csv"), append = TRUE)

message("基因集已保存至: ", GS_DIR)
message("=== 00_setup_genesets.R 运行完毕 ===")
