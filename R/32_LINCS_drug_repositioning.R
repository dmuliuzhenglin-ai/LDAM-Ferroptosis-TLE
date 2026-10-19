#!/usr/bin/env Rscript
# ============================================================================
# Script      : 32_LINCS_drug_repositioning.R
# Purpose     : Drug repositioning against the LDAM-high differential signature via the LINCS L1000CDS2 public API.
# Manuscript  : Methods 2.8 / Results 3.10
# Input       : LDAM-high vs -low DE genes (results/ or recomputed); L1000CDS2 API
# Output      : results/32_* repositioning tables
# Dependencies: Seurat, UCell, httr, jsonlite, tidyverse
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================
# ===========================================================================
# 32_LINCS_drug_repositioning.R  【v2 已实测 API 可用】
# 收官分析 6：基于 LDAM-high 差异签名的药物重定位（L1000CDS2 公共 API）
# 输入：LDAM-high vs LDAM-low 差异基因（优先用 results/ 下已有差异表，
#       找不到则用 FindMarkers 现算），上下调各 top 200
# API： https://amp.pharm.mssm.edu/L1000CDS2/query
#       注意 config 必须带 db-version="latest"，结果直接在响应的 topMeta 里
#       （无需二次 GET）；aggravate=FALSE => score 越负 = 逆转签名越强
# 依赖：Seurat, httr, jsonlite, tidyverse
# ===========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(httr)
  library(jsonlite)
  library(tidyverse)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
setwd(BASE_DIR)
dir.create("results/32_LINCS", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. 差异签名（上下调各 top 200）
# ---------------------------------------------------------------------------
f_de <- list.files("results", pattern = "LDAM.*(high|High).*\\.(csv|txt)$",
                   recursive = TRUE, full.names = TRUE)
f_de <- f_de[!grepl("proportion|score|LINCS", basename(f_de), ignore.case = TRUE)][1]

if (!is.na(f_de)) {
  message("使用已有差异表: ", f_de)
  de <- read_csv(f_de, show_col_types = FALSE)
# [debug] print(colnames(de))
  gcol  <- names(de)[sapply(de, is.character)][1]
  fccol <- names(de)[grepl("log2FC|avg_log2FC|logFC", names(de), ignore.case = TRUE)][1]
  stopifnot(!is.na(gcol), !is.na(fccol))
  up_genes   <- de %>% filter(.data[[fccol]] > 0) %>% arrange(desc(.data[[fccol]])) %>% slice_head(n = 200) %>% pull(gcol)
  down_genes <- de %>% filter(.data[[fccol]] < 0) %>% arrange(.data[[fccol]])        %>% slice_head(n = 200) %>% pull(gcol)
} else {
  message("未找到差异表，用 FindMarkers 现算")
  library(UCell)
  obj <- readRDS("results/06_GSE190452_annotated_full.rds")
  ann_col <- intersect(c("celltype","cell_type","annotation","cluster"), colnames(obj@meta.data))[1]
  mg <- subset(obj, cells = rownames(obj@meta.data)[grepl("micro|MG", obj@meta.data[[ann_col]], ignore.case = TRUE)])
  LDAM_ref <- c("PLIN2","DGAT1","DGAT2","ACAT1","ACSL1","CD36","FABP3","HILPDA","SOAT1","LPIN1")
  mg <- AddModuleScore_UCell(mg, features = list(LDAM_ref = LDAM_ref), ncores = 1)
  mg$LDAM_group <- ifelse(mg$LDAM_ref_UCell > median(mg$LDAM_ref_UCell), "LDAM_high", "LDAM_low")
  Idents(mg) <- mg$LDAM_group
  de <- FindMarkers(mg, ident.1 = "LDAM_high", ident.2 = "LDAM_low",
                    logfc.threshold = 0.1, min.pct = 0.1)
  de$gene <- rownames(de)
  write_csv(de, "results/32_LINCS/32_LDAM_high_vs_low_DE.csv")
  up_genes   <- de %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 200) %>% pull(gene)
  down_genes <- de %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC)        %>% slice_head(n = 200) %>% pull(gene)
}
message("上调签名 ", length(up_genes), " 基因 | 下调签名 ", length(down_genes), " 基因")

# ---------------------------------------------------------------------------
# 2. 提交 L1000CDS2（结果直接在响应里）
# ---------------------------------------------------------------------------
payload <- list(
  data   = list(upGenes = up_genes, dnGenes = down_genes),
  config = list(aggravate = FALSE, searchMethod = "geneSet", share = TRUE,
                combination = FALSE, `db-version` = "latest"),
  meta   = list(title = "LDAM_high_microglia_signature")
)

resp <- POST("https://amp.pharm.mssm.edu/L1000CDS2/query",
             body = toJSON(payload, auto_unbox = TRUE),
             content_type_json(), encode = "raw", timeout(120))
stop_for_status(resp)
job <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE)

top <- job$topMeta
message("返回 ", nrow(top), " 条 L1000 扰动记录")
out <- top %>%
  transmute(perturbagen = pert_desc,
            score       = score,
            cell_line   = cell_id,
            dose        = pert_dose,
            time_h      = pert_time,
            sig_id      = sig_id) %>%
  arrange(score)   # 最负 = 逆转最强
write_csv(out, "results/32_LINCS/32_L1000CDS2_reversal_ranked.csv")

# ---------------------------------------------------------------------------
# 3. 去重汇总：每个化合物取最强逆转记录
# ---------------------------------------------------------------------------
per_cmpd <- out %>%
  filter(!perturbagen %in% c("-666", ""), !is.na(perturbagen)) %>%
  group_by(perturbagen) %>%
  summarise(best_score = min(score), n_records = n(),
            cell_lines = paste(unique(cell_line), collapse = ";"), .groups = "drop") %>%
  arrange(best_score)
write_csv(per_cmpd, "results/32_LINCS/32_L1000CDS2_per_compound.csv")
# [debug] print(head(per_cmpd, 40), n = 40)

message("完成。判读：best_score 越负 = 逆转 LDAM-high 签名越强。")
message("关注：PPAR 激动剂（pioglitazone/troglitazone）、HDAC 抑制剂、")
message("铁死亡相关（ferroptosis 抑制剂少见，但 NRF2 激活剂有意义）、糖皮质激素。")
message("与脚本 30 的 TF 活性结果互相印证：PPARG 活性高 + PPAR 激动剂命中 => 闭环。")
