# ============================================================================
# Script      : FigureS5.R
# Purpose     : Forest plot of Galectin-3 genetic evidence: plasma pQTL + cortical eQTL MR + SMR/HEIDI (Figure S5).
# Manuscript  : Figure S5
# Input       : MR/SMR result tables from 15-17
# Output      : results/FigureS5_MR_forest.png/.pdf
# Dependencies: dplyr, ggplot2, patchwork
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================

## ============================================================
## 18_FigureS5_MR_forest.R
## Galectin-3 遗传证据森林图：血浆 pQTL + 皮层 eQTL MR + SMR/HEIDI
## 输出：results/FigureS5_MR_forest.png / .pdf
## ============================================================
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# install.packages("patchwork")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")

## ---- 配色（与论文一致）----
PURPLE <- "#684765"; PINK <- "#C9828B"; PAPER <- "#F7F3EC"
INK    <- "#2B2430"; GRAY <- "#9B9489"

## ---- 1. 数据（均为脚本15/16/17已验证结果）----
## A: 血浆 pQTL（OR 尺度）；focal/all 为阴性（仅 P 值，无 CI 展示）
dfA <- data.frame(
  outcome = factor(c("Focal epilepsy + HS\n(ieu-b-14)",
                     "Focal epilepsy, all\n(ieu-b-10)",
                     "All epilepsy\n(ieu-b-8)"),
                   levels = c("All epilepsy\n(ieu-b-8)",
                              "Focal epilepsy, all\n(ieu-b-10)",
                              "Focal epilepsy + HS\n(ieu-b-14)")),
  est  = c(1.0071, NA, NA),
  lo   = c(1.0021, NA, NA),
  hi   = c(1.0122, NA, NA),
  plab = c("P = 5.6e-03", "P = 0.31 (n.s.)", "P = 0.36 (n.s.)"),
  sig  = c(TRUE, FALSE, FALSE)
)

## B: 皮层 eQTL MR（beta 尺度，脚本16结果）
dfB <- data.frame(
  outcome = factor(c("Focal epilepsy + HS\n(ieu-b-14)",
                     "Focal epilepsy, all\n(ieu-b-10)",
                     "All epilepsy\n(ieu-b-8)"),
                   levels = levels(dfA$outcome)),
  est  = c(-0.01582, -0.02514, -0.03618),
  se   = c(0.004382, 0.025902, 0.020688),
  plab = c("P = 3.1e-04", "P = 0.33", "P = 0.080"),
  sig  = c(TRUE, FALSE, FALSE)
) %>% mutate(lo = est - 1.96 * se, hi = est + 1.96 * se)

## ---- 2. 面板 A：血浆 pQTL ----
pA <- ggplot(dfA, aes(y = outcome)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = GRAY, linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, color = sig), height = 0, linewidth = 1.1, na.rm = TRUE) +
  geom_point(aes(x = est, color = sig), shape = 18, size = 4.5, na.rm = TRUE) +
  geom_point(data = subset(dfA, !sig), x = 1, shape = 21, size = 3.5,
             fill = "white", color = GRAY, stroke = 1) +
  geom_text(aes(x = 1.0148, label = plab, color = sig, fontface = ifelse(sig, "bold", "plain")),
            hjust = 0, size = 3.4, show.legend = FALSE) +
  annotate("text", x = 0.9973, y = 3.4, label = "subtype-specific",
           color = PINK, fontface = "italic", size = 3.2, hjust = 0) +
  scale_color_manual(values = c("TRUE" = PURPLE, "FALSE" = GRAY), guide = "none") +
  scale_x_continuous(limits = c(0.997, 1.0225)) +
  coord_cartesian(ylim = c(0.6, 3.6), clip = "off") +
  labs(title = "A  Plasma Galectin-3 pQTL → epilepsy",
       x = "OR per SD plasma Galectin-3 (IVW)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, color = INK),
        axis.text = element_text(color = INK, size = 10),
        axis.title.x = element_text(color = INK, size = 10),
        panel.grid = element_blank(),
        axis.line = element_line(color = GRAY),
        plot.background = element_rect(fill = PAPER, color = NA),
        panel.background = element_rect(fill = PAPER, color = NA))

## ---- 3. 面板 B：皮层 eQTL MR ----
pB <- ggplot(dfB, aes(y = outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = GRAY, linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, color = sig), height = 0, linewidth = 1.1) +
  geom_point(aes(x = est, color = sig, shape = sig), size = 4.5) +
  geom_text(aes(x = 0.033, label = plab, color = sig, fontface = ifelse(sig, "bold", "plain")),
            hjust = 0, size = 3.4, show.legend = FALSE) +
  annotate("text", x = -0.096, y = 3.4, label = "protective direction — mirror of plasma",
           color = PINK, fontface = "italic", size = 3.2, hjust = 0) +
  scale_color_manual(values = c("TRUE" = PURPLE, "FALSE" = GRAY), guide = "none") +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
  scale_x_continuous(limits = c(-0.098, 0.062)) +
  coord_cartesian(ylim = c(0.6, 3.6), clip = "off") +
  labs(title = "B  Cortex LGALS3 eQTL → epilepsy (MetaBrain)",
       x = "β per SD cortical LGALS3 expression (IVW)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, color = INK),
        axis.text = element_text(color = INK, size = 10),
        axis.title.x = element_text(color = INK, size = 10),
        panel.grid = element_blank(),
        axis.line = element_line(color = GRAY),
        plot.background = element_rect(fill = PAPER, color = NA),
        panel.background = element_rect(fill = PAPER, color = NA))

## ---- 4. 拼接 + SMR 注释条 ----
library(patchwork)

smr_note <- ggplot() +
  annotate("label", x = 0.5, y = 0.5,
           label = paste0("SMR (cis, top SNP rs8021940): b = -0.011, P(SMR) = 0.068   |   ",
                          "multi-SNP SMR: b = -0.010, P = 0.087   |   ",
                          "HEIDI P = 0.51 / 0.36 (no LD confounding)"),
           size = 3.2, color = INK, fill = "white", label.r = unit(0.15, "lines"),
           label.padding = unit(0.4, "lines"), linewidth = 0.4) +
  xlim(0, 1) + ylim(0, 1) + theme_void() +
  theme(plot.background = element_rect(fill = PAPER, color = NA))
## 注释条边框用粉色
smr_note$layers[[1]]$aes_params$colour <- PINK

figS5 <- (pA | pB) / smr_note + plot_layout(heights = c(10, 1.2))
figS5 <- figS5 & theme(plot.background = element_rect(fill = PAPER, color = NA))

ggsave(file.path(OUT_DIR, "FigureS5_MR_forest.png"), figS5,
       width = 13, height = 5.4, dpi = 300, bg = PAPER)
ggsave(file.path(OUT_DIR, "FigureS5_MR_forest.pdf"), figS5,
       width = 13, height = 5.4, bg = PAPER)

cat("已保存:\n", file.path(OUT_DIR, "FigureS5_MR_forest.png"), "\n",
    file.path(OUT_DIR, "FigureS5_MR_forest.pdf"), "\n")
