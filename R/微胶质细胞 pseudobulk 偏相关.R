# ============================================================================
# Script      : 微胶质细胞 pseudobulk 偏相关.R
# Purpose     : Sample-level pseudobulk partial correlation: LGALS3 x ferroptosis controlling LDAM (with the reverse control).
# Manuscript  : Methods 2.6 / Results 3.8 (Figure S7)
# Input       : GSE190452 microglia pseudobulk (results/06)
# Output      : results/20_pseudobulk_partial_cor.csv + Figure S7
# Dependencies: Seurat, dplyr, ggplot2, ppcor
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================

## ============================================================
## 20_pseudobulk_partial_cor.R（整段一次性运行）
## 样本级 pseudobulk：LGALS3 × 铁死亡 偏相关（控制 LDAM 分数）
## 对照：LGALS3 × LDAM（控制铁死亡）应为阴性 → 双向证明特异性
## 产出：results/20_pseudobulk_partial_cor.csv + FigureS7
## ============================================================
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# install.packages("ppcor")
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(ppcor)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")
PURPLE <- "#684765"; PINK <- "#C9828B"; PAPER <- "#F7F3EC"; INK <- "#2B2430"; GRAY <- "#9B9489"

## ---- 1. 指定两个微胶质对象 ----
RDS_FILES <- c(
  GSE190452 = file.path(BASE_DIR, "results", "06_GSE190452_microglia.rds"),
  GSE201048 = file.path(BASE_DIR, "results", "07_GSE201048_microglia.rds")
)

FERRO <- c("SLC7A11","GPX4","GCLM","GSTP1","HSPB1","ACSL4","CP","FTL1","FTH1",
           "ALOX15","PTGS2","NQO1","HMOX1","SAT1","TFRC")
LDAM  <- c("APOE","CH25H","ACSL1","C1QA","C1QB","C1QC","TYROBP","CD9","PLIN2",
           "SPP1","CD44","LPL","CST7","TREM2","AXL")

run_one <- function(rds, tag, sample_col) {
  cat("\n==========", tag, "==========\n")
  mg <- readRDS(rds)
  meta <- mg@meta.data
  cat("细胞数:", ncol(mg), "\n")
  stopifnot(sample_col %in% colnames(meta))
  cat("样本列:", sample_col, " | 样本数:", length(unique(meta[[sample_col]])), "\n")

  ## 优先复用已有 UCell 分数；没有再算 AddModuleScore
  if (all(c("LDAM_UCell", "Ferroptosis_UCell") %in% colnames(meta))) {
    ldam_v <- meta$LDAM_UCell; ferro_v <- meta$Ferroptosis_UCell
    cat("使用已有 LDAM_UCell / Ferroptosis_UCell\n")
  } else {
    mg <- AddModuleScore(mg, features = list(intersect(FERRO, rownames(mg))), name = "FERRO", seed = 1)
    mg <- AddModuleScore(mg, features = list(intersect(LDAM,  rownames(mg))), name = "LDAM",  seed = 1)
    ldam_v <- mg$LDAM1; ferro_v <- mg$FERRO1
    cat("使用 AddModuleScore 重新计算\n")
  }

  det <- mean(FetchData(mg, "LGALS3")[,1] > 0)
  cat("LGALS3 检出率:", sprintf("%.2f%%", 100*det), "\n")

  pb <- data.frame(sample = meta[[sample_col]],
                   LGALS3 = FetchData(mg, "LGALS3")[,1],
                   FERRO = ferro_v, LDAM = ldam_v) %>%
    group_by(sample) %>%
    summarise(LGALS3 = mean(LGALS3), FERRO = mean(FERRO), LDAM = mean(LDAM),
              n_cells = n(), .groups = "drop")
  cat("pseudobulk 样本:\n"); print(as.data.frame(pb))

  if (nrow(pb) < 4) { cat("样本数 <4，跳过偏相关\n"); return(list(pb = pb, res = NULL)) }
  res <- data.frame(
    dataset = tag,
    test = c("LGALS3~FERRO | LDAM", "LGALS3~LDAM | FERRO", "FERRO~LDAM | LGALS3"),
    r = c(pcor.test(pb$LGALS3, pb$FERRO, pb$LDAM)$estimate,
          pcor.test(pb$LGALS3, pb$LDAM,  pb$FERRO)$estimate,
          pcor.test(pb$FERRO,  pb$LDAM,  pb$LGALS3)$estimate),
    P = c(pcor.test(pb$LGALS3, pb$FERRO, pb$LDAM)$p.value,
          pcor.test(pb$LGALS3, pb$LDAM,  pb$FERRO)$p.value,
          pcor.test(pb$FERRO,  pb$LDAM,  pb$LGALS3)$p.value),
    n_samples = nrow(pb), LGALS3_detection = det)
  print(res)
  list(pb = pb, res = res)
}

out1 <- run_one(RDS_FILES["GSE190452"], "GSE190452", sample_col = "orig.ident")
out2 <- run_one(RDS_FILES["GSE201048"], "GSE201048", sample_col = "sampleName")
res_all <- bind_rows(out1$res, out2$res)
if (nrow(res_all)) write.csv(res_all, file.path(OUT_DIR, "20_pseudobulk_partial_cor.csv"), row.names = FALSE)
# [deprecated] ---- orphaned remnant of an older function version (invalid in source); kept commented for provenance ----
# [deprecated]   ## 模块分数
# [deprecated]   mg <- AddModuleScore(mg, features = list(intersect(FERRO, rownames(mg))), name = "FERRO", seed = 1)
# [deprecated]   mg <- AddModuleScore(mg, features = list(intersect(LDAM,  rownames(mg))), name = "LDAM",  seed = 1)
# [deprecated]   ## ---- Figure S7：GSE190452 样本级三联图 ----
# [deprecated]   pb <- out1$pb %>%
# [deprecated]     mutate(group = ifelse(grepl("TLE", sample), "TLE", "Normal"))

# [deprecated]   p <- ggplot(pb, aes(LDAM, FERRO)) +
# [deprecated]     geom_point(aes(size = LGALS3, color = LGALS3), alpha = 0.9) +
# [deprecated]     geom_smooth(method = "lm", se = FALSE, color = GRAY, linetype = "dashed", linewidth = 0.6) +
# [deprecated]     ggrepel::geom_text_repel(aes(label = gsub("GSM[0-9]+_|_matrix.mtx.gz", "", sample)),
# [deprecated]                              size = 3.2, color = INK) +
# [deprecated]     scale_color_gradient(low = GRAY, high = PURPLE, name = "LGALS3\n(mean expr)") +
# [deprecated]     scale_size(range = c(3, 9), name = "LGALS3\n(mean expr)") +
# [deprecated]     labs(title = "Microglia sample-level pseudobulk (GSE190452, n = 8)",
# [deprecated]          subtitle = "partial r(LGALS3, FERRO | LDAM) = 0.821, P = 0.024;  partial r(LGALS3, LDAM | FERRO) = -0.049, P = 0.92",
# [deprecated]          x = "LDAM score (UCell)", y = "Ferroptosis score (UCell)") +
# [deprecated]     theme_minimal(base_size = 11) +
# [deprecated]     theme(plot.title = element_text(face = "bold", size = 12, color = INK),
# [deprecated]           plot.subtitle = element_text(size = 9, color = PURPLE),
# [deprecated]           text = element_text(color = INK),
# [deprecated]           plot.background = element_rect(fill = PAPER, color = NA),
# [deprecated]           panel.background = element_rect(fill = PAPER, color = NA))
# [deprecated]   ggsave(file.path(OUT_DIR, "FigureS7_pseudobulk_partial.png"), p,
# [deprecated]          width = 7.5, height = 5.8, dpi = 300, bg = PAPER)
# [deprecated]   ggsave(file.path(OUT_DIR, "FigureS7_pseudobulk_partial.pdf"), p,
# [deprecated]          width = 7.5, height = 5.8, bg = PAPER)
# [deprecated]   print(p)
# [deprecated]   ## LGALS3 检出率（太低则该数据集 mRNA 层面不可做——这正是论文"蛋白层面诱导"的论点）
# [deprecated]   det <- mean(FetchData(mg, "LGALS3")[,1] > 0)
# [deprecated]   cat("LGALS3 检出率:", sprintf("%.2f%%", 100*det), "\n")

# [deprecated]   pb <- data.frame(
# [deprecated]     sample = meta[[sample_col]],
# [deprecated]     LGALS3 = FetchData(mg, "LGALS3")[,1],
# [deprecated]     FERRO  = mg$FERRO1,
# [deprecated]     LDAM   = mg$LDAM1
# [deprecated]   ) %>% group_by(sample) %>%
# [deprecated]     summarise(LGALS3 = mean(LGALS3), FERRO = mean(FERRO), LDAM = mean(LDAM),
# [deprecated]               n_cells = n(), .groups = "drop")
# [deprecated]   cat("pseudobulk 样本:\n"); print(as.data.frame(pb))

# [deprecated]   res <- data.frame(
# [deprecated]     dataset = tag,
# [deprecated]     test = c("LGALS3~FERRO | LDAM", "LGALS3~LDAM | FERRO", "FERRO~LDAM | LGALS3"),
# [deprecated]     r = c(pcor.test(pb$LGALS3, pb$FERRO, pb$LDAM)$estimate,
# [deprecated]           pcor.test(pb$LGALS3, pb$LDAM,  pb$FERRO)$estimate,
# [deprecated]           pcor.test(pb$FERRO,  pb$LDAM,  pb$LGALS3)$estimate),
# [deprecated]     P = c(pcor.test(pb$LGALS3, pb$FERRO, pb$LDAM)$p.value,
# [deprecated]           pcor.test(pb$LGALS3, pb$LDAM,  pb$FERRO)$p.value,
# [deprecated]           pcor.test(pb$FERRO,  pb$LDAM,  pb$LGALS3)$p.value),
# [deprecated]     n_samples = nrow(pb),
# [deprecated]     LGALS3_detection = det)
# [deprecated]   print(res)
# [deprecated]   list(pb = pb, res = res)
# [deprecated] }

out1 <- run_one(RDS_FILES["GSE190452"], "GSE190452")
out2 <- run_one(RDS_FILES["GSE201048"], "GSE201048")
res_all <- bind_rows(out1$res, out2$res)
write.csv(res_all, file.path(OUT_DIR, "20_pseudobulk_partial_cor.csv"), row.names = FALSE)
cat("\n已保存:", file.path(OUT_DIR, "20_pseudobulk_partial_cor.csv"), "\n")
## ---- 6. Figure S7 ----
p <- ggplot(pb, aes(LDAM, FERRO, size = LGALS3, color = LGALS3)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = GRAY, high = PURPLE) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3, color = INK) +
  labs(title = "Sample-level pseudobulk: microglia",
       subtitle = sprintf("partial r(LGALS3,FERRO|LDAM) = %.3f, P = %.3g",
                          res$r[1], res$P[1]),
       x = "LDAM signature score", y = "Ferroptosis score", color = "LGALS3", size = "LGALS3") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", color = INK),
        text = element_text(color = INK),
        plot.background = element_rect(fill = PAPER, color = NA),
        panel.background = element_rect(fill = PAPER, color = NA))
ggsave(file.path(OUT_DIR, "FigureS7_pseudobulk_partial.png"), p,
       width = 7, height = 5.5, dpi = 300, bg = PAPER)
cat("\n已保存 results/20_pseudobulk_partial_cor.csv 与 FigureS7\n")
