# ============================================================================
# Script      : 扰动拯救分析.R
# Purpose     : Perturbation-rescue re-analysis (JCI 2023, PMC10786694): Lgals3 KO in Tau22 mouse hippocampus; TD139 in human iMGL.
# Manuscript  : Methods 2.7 / Results 3.9 (Figure S6)
# Input       : JCI2023 *.xlsx tables
# Output      : results/19_rescue_summary.csv + Figure S6
# Dependencies: dplyr, ggplot2, patchwork, readxl
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 19_Lgals3_perturbation_rescue.R（最终版，整段一次性运行）
## JCI 2023 (PMC10786694) 公开扰动数据重分析：
##   A. 小鼠 Tau22 海马：Lgals3 KO 对 Tau22 转录组的拯救
##   B. 人源 iMGL：TD139（Gal-3 抑制剂）对 pTau 转录组的拯救
## 产出：results/19_rescue_summary.csv + FigureS6_rescue.png/pdf
## ============================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(ggplot2); library(patchwork)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")

F_MOUSE <- file.path(BASE_DIR, "data", "JCI2023_Tau22_Lgals3KO_差异基因表.xlsx")
F_HUMAN <- file.path(BASE_DIR, "data", "JCI2023_iMGL_TD139_差异基因表.xlsx")

PURPLE <- "#684765"; PINK <- "#C9828B"; PAPER <- "#F7F3EC"; INK <- "#2B2430"; GRAY <- "#9B9489"

## ---- 1. 小鼠：Tau22/Lgals3-WT vs WT（疾病效应）+ Tau22/Lgals3-KO vs Tau22/WT（拯救效应）----
wt <- read_excel(F_MOUSE, sheet = "6B; TG3WT_vs_WG3WT", skip = 1)
ko <- read_excel(F_MOUSE, sheet = "6B; TG3KO_vs_TG3WT", skip = 1)
colnames(wt) <- c("gene_id","gene","ncbi","log2fc_Tau22","p_Tau22")
colnames(ko) <- c("gene_id","gene","ncbi","log2fc_KO","p_KO")

m <- inner_join(wt, ko, by = c("gene_id","gene"))
m$reversed <- sign(m$log2fc_Tau22) * sign(m$log2fc_KO) < 0

rescue_mouse <- binom.test(sum(m$reversed), nrow(m), p = 0.5)
cor_mouse    <- cor.test(m$log2fc_Tau22, m$log2fc_KO, method = "spearman")
cat(sprintf("小鼠: %d/%d (%.1f%%) 被KO逆转, binom P = %.2e, Spearman rho = %.3f (P = %.2e)\n",
            sum(m$reversed), nrow(m), 100*mean(m$reversed),
            rescue_mouse$p.value, cor_mouse$estimate, cor_mouse$p.value))

## ---- 2. 人源 iMGL：pTau vs Con（损伤）+ TD139 vs pTau（拯救）----
ptau <- read_excel(F_HUMAN, sheet = "2B; pT vs Con Log2 0.55", skip = 1)
td   <- read_excel(F_HUMAN, sheet = "2B; pTauTD139 vs pTau Log2 0.55", skip = 1)
colnames(ptau) <- c("gene_id","gene","ncbi","log2fc_pTau","p_pTau")
colnames(td)   <- c("gene_id","gene","ncbi","log2fc_TD139","p_TD139")

h <- inner_join(ptau, td, by = c("gene_id","gene"))
h$reversed <- sign(h$log2fc_pTau) * sign(h$log2fc_TD139) < 0

rescue_human <- binom.test(sum(h$reversed), nrow(h), p = 0.5)
cor_human    <- cor.test(h$log2fc_pTau, h$log2fc_TD139, method = "spearman")
cat(sprintf("人iMGL: %d/%d (%.1f%%) 被TD139逆转, binom P = %.2e, Spearman rho = %.3f (P = %.2e)\n",
            sum(h$reversed), nrow(h), 100*mean(h$reversed),
            rescue_human$p.value, cor_human$estimate, cor_human$p.value))

## ---- 3. 关键通路基因命中 ----
FERRO <- c("SLC7A11","GPX4","GCLM","GSTP1","HSPB1","ACSL4","CP","FTL1","FTH1",
           "ALOX15","PTGS2","NQO1","HMOX1","SAT1","TFRC")
LDAM  <- c("APOE","CH25H","ACSL1","C1QA","C1QB","C1QC","TYROBP","CD9","PLIN2",
           "SPP1","CD44","LPL","CST7","TREM2","AXL")
GLYC  <- c("MGAT1","MGAT2","MGAT4A","MGAT4B","MGAT5","ST3GAL1","ST3GAL4",
           "ST6GAL1","B4GALT1","B4GALT5","OGT","OGA","HEXB","GUSB","NEU1")

to_mouse <- function(v) paste0(toupper(substr(v,1,1)), tolower(substr(v,2,nchar(v))))

hits_ko <- ko %>% filter(gene %in% c(to_mouse(FERRO), to_mouse(LDAM), to_mouse(GLYC))) %>%
  mutate(set = case_when(gene %in% to_mouse(FERRO) ~ "ferroptosis",
                         gene %in% to_mouse(LDAM)  ~ "LDAM", TRUE ~ "glycosylation"))
hits_td <- td %>% filter(gene %in% c(FERRO, LDAM, GLYC)) %>%
  mutate(set = case_when(gene %in% FERRO ~ "ferroptosis",
                         gene %in% LDAM  ~ "LDAM", TRUE ~ "glycosylation"))
cat("\n[KO 表关键基因]\n");   print(as.data.frame(hits_ko))
cat("\n[TD139 表关键基因]\n"); print(as.data.frame(hits_td))

## ---- 4. 汇总表 ----
summary_df <- data.frame(
  model   = c("Mouse Tau22 hippocampus (Lgals3-KO)", "Human iMGL + pTau (TD139)"),
  n_shared_DEG = c(nrow(m), nrow(h)),
  pct_reversed = c(mean(m$reversed), mean(h$reversed)),
  binom_P = c(rescue_mouse$p.value, rescue_human$p.value),
  spearman_rho = c(unname(cor_mouse$estimate), unname(cor_human$estimate)),
  cor_P = c(cor_mouse$p.value, cor_human$p.value))
write.csv(summary_df, file.path(OUT_DIR, "19_rescue_summary.csv"), row.names = FALSE)

## ---- 5. Figure S6：双面板拯救散点图 ----
mk_panel <- function(df, x, y, rev_col, xname, yname, title, key_genes) {
  df$key <- df$gene %in% key_genes
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_hline(yintercept = 0, color = GRAY, linewidth = 0.4) +
    geom_vline(xintercept = 0, color = GRAY, linewidth = 0.4) +
    geom_point(aes(color = .data[[rev_col]]), alpha = 0.55, size = 1.4) +
    geom_point(data = subset(df, key), color = INK, size = 2.5) +
    ggrepel::geom_text_repel(data = subset(df, key), aes(label = gene),
                             size = 3.2, color = INK, max.overlaps = 20) +
    scale_color_manual(values = c("TRUE" = PINK, "FALSE" = "#B9B2A8"),
                       labels = c("concordant", "reversed"), name = NULL) +
    labs(title = title, x = xname, y = yname) +
    coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11.5, color = INK),
          legend.position = "top", text = element_text(color = INK),
          plot.background = element_rect(fill = PAPER, color = NA),
          panel.background = element_rect(fill = PAPER, color = NA))
}

pA <- mk_panel(m, "log2fc_Tau22", "log2fc_KO", "reversed",
               "log2FC (Tau22 vs WT)", "log2FC (Lgals3-KO vs Tau22)",
               "A  Mouse Tau22 hippocampus: Lgals3 knockout rescue",
               to_mouse(c("CD44","HSPB1","NQO1")))
pB <- mk_panel(h, "log2fc_pTau", "log2fc_TD139", "reversed",
               "log2FC (pTau vs Con)", "log2FC (TD139 + pTau vs pTau)",
               "B  Human iMGL: Gal-3 inhibitor (TD139) rescue",
               c("TFRC","CP"))

figS6 <- pA + pB + plot_annotation(
  caption = sprintf("Reversed: mouse %d/%d (%.0f%%), P = %.1e; human iMGL %d/%d (%.0f%%), P = %.1e",
                    sum(m$reversed), nrow(m), 100*mean(m$reversed), rescue_mouse$p.value,
                    sum(h$reversed), nrow(h), 100*mean(h$reversed), rescue_human$p.value),
  theme = theme(plot.caption = element_text(color = INK, size = 9.5),
                plot.background = element_rect(fill = PAPER, color = NA)))

ggsave(file.path(OUT_DIR, "FigureS6_rescue.png"), figS6, width = 12, height = 5.6, dpi = 300, bg = PAPER)
ggsave(file.path(OUT_DIR, "FigureS6_rescue.pdf"), figS6, width = 12, height = 5.6, bg = PAPER)
cat("\n已保存 FigureS6 与 19_rescue_summary.csv 到 results/\n")
