# ============================================================================
# Script      : Gal-3 受体轴分析.R
# Purpose     : Gal-3 lectin-glycan-receptor axis: receptor expression in LDAM-high microglia, co-expression, CellChat LGALS3->receptor signalling, perturbation response.
# Manuscript  : Methods 2.6 / Results 3.8 (Figure S9)
# Input       : results/06 rds; CellChat objects; JCI2023 perturbation tables
# Output      : results/22_*.csv + Figure S9
# Dependencies: Seurat, dplyr, ggplot2, ggrepel, patchwork, readxl
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 22_Gal3_receptor_axis.R（整段一次性运行）
## Gal-3 凝集素-糖链-受体轴的四层证据：
##   A. Gal-3 受体在 LDAM-high vs LDAM-low 微胶质中的表达
##   B. 受体/糖链合成酶与 LGALS3 的共表达
##   C. CellChat 自定义对象中 LGALS3→受体通讯强度（TLE vs Control）
##   D. 受体基因在 Lgals3-KO / TD139 扰动中的变化（复用19的数据）
## 产出：results/22_*.csv + FigureS9
## ============================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(readxl); library(ggrepel)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")
PURPLE <- "#684765"; PINK <- "#C9828B"; PAPER <- "#F7F3EC"; INK <- "#2B2430"; GRAY <- "#9B9489"

## Gal-3 已知糖链受体 + poly-LacNAc 合成酶（Gal-3 配体糖链）
RECEPTORS <- c("TREM2","CD44","ITGAM","ITGB2","ITGB1","ITGAX","TLR4","PTPRC","MERTK","CD98"="SLC3A2")
RECEPTORS <- c("TREM2","CD44","ITGAM","ITGB2","ITGB1","ITGAX","TLR4","PTPRC","MERTK","SLC3A2")
GLYCAN_SYN <- c("MGAT5","MGAT5B","B4GALT1","B4GALT2","B4GALT5","ST3GAL1","ST3GAL4")  # poly-LacNAc / N-聚糖分支

## ---- A+B 无函数版：两个数据集循环 ----
resA_all <- list(); resB_all <- list()

for (tag in c("GSE190452", "GSE201048")) {
  rds <- file.path(BASE_DIR, "results",
                   ifelse(tag == "GSE190452", "06_GSE190452_microglia.rds", "07_GSE201048_microglia.rds"))
  cat("\n==========", tag, "==========\n")
  mg <- readRDS(rds)
  if (inherits(mg[["RNA"]], "Assay5")) mg <- JoinLayers(mg)   # 仅 v5 多 layer 对象需要

  if ("LDAM_level" %in% colnames(meta)) {
    high_lab <- grep("high", unique(as.character(meta$LDAM_level)), ignore.case = TRUE, value = TRUE)[1]
    is_high <- as.character(meta$LDAM_level) == high_lab
  } else {
    is_high <- meta$LDAM_UCell > median(meta$LDAM_UCell)
  }
  cat("LDAM-high:", sum(is_high), " LDAM-low:", sum(!is_high), "\n")

  genes_use <- intersect(c(RECEPTORS, GLYCAN_SYN, "LGALS3"), rownames(mg))
  cat("genes_use:", length(genes_use), "\n")
  expr <- as.matrix(GetAssayData(mg, layer = "data")[genes_use, , drop = FALSE])

  resA <- lapply(genes_use, function(g) {
    x <- expr[g, is_high]; y <- expr[g, !is_high]
    p <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA)
    data.frame(dataset = tag, gene = g,
               type = ifelse(g %in% RECEPTORS, "receptor",
                             ifelse(g %in% GLYCAN_SYN, "glycan_synth", "lectin")),
               mean_high = mean(x), mean_low = mean(y),
               log2FC = log2((mean(x) + 1e-6) / (mean(y) + 1e-6)),
               det_high = mean(x > 0), det_low = mean(y > 0), P = p)
  }) %>% bind_rows()
  resA$P_adj <- p.adjust(resA$P, "BH")
  print(as.data.frame(resA[order(resA$P), c("gene","type","log2FC","P","P_adj")]))
  resA_all[[tag]] <- resA

  if ("LGALS3" %in% genes_use && sum(expr["LGALS3", ] > 0) > 20) {
    lg <- expr["LGALS3", is_high]
    resB <- lapply(setdiff(genes_use, "LGALS3"), function(g) {
      ct <- suppressWarnings(cor.test(lg, expr[g, is_high], method = "spearman"))
      data.frame(dataset = tag, gene = g,
                 type = ifelse(g %in% RECEPTORS, "receptor", "glycan_synth"),
                 rho = ct$estimate, P = ct$p.value)
    }) %>% bind_rows()
    resB$P_adj <- p.adjust(resB$P, "BH")
    cat("\n[LDAM-high 内与 LGALS3 共表达 Top10]:\n")
    print(as.data.frame(resB[order(resB$P), ])[1:10, ])
    resB_all[[tag]] <- resB
  } else {
    cat("LGALS3 检出过低，共表达跳过\n")
  }
}

resA_all <- bind_rows(resA_all)
resB_all <- bind_rows(resB_all)
write.csv(resA_all, file.path(OUT_DIR, "22_receptor_LDAM_highlow.csv"), row.names = FALSE)
write.csv(resB_all, file.path(OUT_DIR, "22_receptor_LGALS3_coexpr.csv"), row.names = FALSE)
## ============================================================
## C. CellChat 自定义对象中的 LGALS3 通讯
## ============================================================
cat("\n========== CellChat LGALS3 通讯 ==========\n")
cc_file <- file.path(BASE_DIR, "results", "14_CellChat_gal_custom.rds")
if (file.exists(cc_file)) {
  cc <- readRDS(cc_file)
# [debug] print(names(cc@net))
  ## 通讯概率表（不同版本字段不同，做防御式提取）
  prob <- cc@net$prob
  if (!is.null(prob)) {
    lr_names <- dimnames(prob)[[3]]
    gal_idx <- grep("LGALS3|Gal", lr_names, ignore.case = TRUE)
    cat("含 LGALS3 的 LR 对:", lr_names[gal_idx], "\n")
    if (length(gal_idx)) {
      comm <- lapply(gal_idx, function(k) {
        m <- prob[,,k]
        data.frame(LR = lr_names[k],
                   sender = rep(rownames(m), ncol(m)),
                   receiver = rep(colnames(m), each = nrow(m)),
                   prob = as.vector(m))
      }) %>% bind_rows() %>% filter(prob > 0) %>% arrange(-prob)
      print(comm)
      write.csv(comm, file.path(OUT_DIR, "22_cellchat_LGALS3_comm.csv"), row.names = FALSE)
    }
  }
  ## 若对象里存了分组比较（TLE vs Control），提取 netP 通路强度
  if (!is.null(cc@netP$pathways)) {
    gal_pw <- grep("LGALS3|GALECTIN|Gal", cc@netP$pathways, value = TRUE, ignore.case = TRUE)
    cat("含 Gal 的通路:", gal_pw, "\n")
  }
} else cat("未找到", cc_file, "，跳过\n")

## ============================================================
## D. 扰动数据中的受体变化（复用19的文件）
## ============================================================
cat("\n========== 扰动数据中的受体 ==========\n")
F_MOUSE <- file.path(BASE_DIR, "data", "JCI2023_Tau22_Lgals3KO_差异基因表.xlsx")
F_HUMAN <- file.path(BASE_DIR, "data", "JCI2023_iMGL_TD139_差异基因表.xlsx")
to_mouse <- function(v) paste0(toupper(substr(v,1,1)), tolower(substr(v,2,nchar(v))))

ko <- read_excel(F_MOUSE, sheet = "6B; TG3KO_vs_TG3WT", skip = 1)
colnames(ko) <- c("gene_id","gene","ncbi","log2fc","p")
td <- read_excel(F_HUMAN, sheet = "2B; pTauTD139 vs pTau Log2 0.55", skip = 1)
colnames(td) <- c("gene_id","gene","ncbi","log2fc","p")

d_ko <- ko %>% filter(gene %in% c(to_mouse(RECEPTORS), to_mouse(GLYCAN_SYN))) %>%
  mutate(dataset = "Mouse Lgals3-KO", direction = ifelse(log2fc < 0, "down", "up"))
d_td <- td %>% filter(gene %in% c(RECEPTORS, GLYCAN_SYN)) %>%
  mutate(dataset = "Human iMGL TD139", direction = ifelse(log2fc < 0, "down", "up"))
resD <- bind_rows(d_ko, d_td)
print(as.data.frame(resD))
write.csv(resD, file.path(OUT_DIR, "22_receptor_perturbation.csv"), row.names = FALSE)

## ============================================================
## Figure S9：受体差异热图（A）+ 扰动变化（B）
## ============================================================
pA <- ggplot(resA_all, aes(dataset, gene, fill = log2FC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(P_adj < 0.05, "*", "")), color = INK, size = 6) +
  scale_fill_gradient2(low = PINK, mid = PAPER, high = PURPLE, midpoint = 0,
                       name = "log2FC\n(LDAM-high/low)") +
  labs(title = "A  Gal-3 receptors & glycan synthases in LDAM-high microglia",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(text = element_text(color = INK), plot.title = element_text(face = "bold", size = 11.5),
        plot.background = element_rect(fill = PAPER, color = NA),
        panel.background = element_rect(fill = PAPER, color = NA),
        axis.text.x = element_text(angle = 0))

if (nrow(resD)) {
  pB <- ggplot(resD, aes(reorder(gene, log2fc), log2fc, fill = direction)) +
    geom_col(width = 0.7) + coord_flip() +
    facet_wrap(~dataset, scales = "free_y", ncol = 1) +
    geom_hline(yintercept = 0, color = GRAY, linewidth = 0.5) +
    scale_fill_manual(values = c(down = PINK, up = PURPLE), name = NULL) +
    labs(title = "B  Receptor/glycan genes upon Gal-3 perturbation",
         x = NULL, y = "log2FC") +
    theme_minimal(base_size = 11) +
    theme(text = element_text(color = INK), plot.title = element_text(face = "bold", size = 11.5),
          plot.background = element_rect(fill = PAPER, color = NA),
          panel.background = element_rect(fill = PAPER, color = NA))
} else {
  pB <- ggplot() + annotate("text", 0.5, 0.5, label = "no receptor hits") + theme_void()
}

figS9 <- (pA | pB) + plot_layout(widths = c(1.6, 1))
ggsave(file.path(OUT_DIR, "FigureS9_Gal3_receptor_axis.png"), figS9,
       width = 12.5, height = 6, dpi = 300, bg = PAPER)
ggsave(file.path(OUT_DIR, "FigureS9_Gal3_receptor_axis.pdf"), figS9,
       width = 12.5, height = 6, bg = PAPER)
cat("\n已保存 results/22_*.csv 与 FigureS9\n")
