# ============================================================================
# Script      : NicheNet 配体-靶标分析.R
# Purpose     : NicheNet ligand-target prediction: microglial LGALS3 -> neurons; ferroptosis-program enrichment of predicted targets.
# Manuscript  : Methods 2.6 related supplementary (Figure S8)
# Input       : GSE190452 Seurat object; nichenet models
# Output      : results/21_nichenet_*.csv + Figure S8
# Dependencies: Seurat, dplyr, ggplot2, ggrepel, nichenetr, patchwork
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================

## ============================================================
## 21_NicheNet_LGALS3.R（整段一次性运行）
## 小胶质细胞 LGALS3(配体) → 神经元(接收) 下游靶基因预测
## 检验：LGALS3 的预测靶基因是否富集铁死亡程序
## 产出：results/21_nichenet_*.csv + FigureS8
## 依赖：nichenetr（首次需安装，见下）
## ============================================================

## 首次安装（已装可跳过）：
##
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# install.packages("nichenetr")
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# remotes::install_github("saeyslab/nichenetr")

suppressPackageStartupMessages({
  library(Seurat); library(nichenetr); library(dplyr)
  library(ggplot2); library(ggrepel)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")
PURPLE <- "#684765"; PINK <- "#C9828B"; PAPER <- "#F7F3EC"; INK <- "#2B2430"; GRAY <- "#9B9489"

FERRO <- c("SLC7A11","GPX4","GCLM","GSTP1","HSPB1","ACSL4","CP","FTL1","FTH1",
           "ALOX15","PTGS2","NQO1","HMOX1","SAT1","TFRC")

## ---- 1. 下载 NicheNet 模型文件（人源，只需下一次）----
MODEL_DIR <- file.path(BASE_DIR, "nichenet_models"); dir.create(MODEL_DIR, showWarnings = FALSE)
lt_file <- file.path(MODEL_DIR, "ligand_target_matrix.rds")          # 改这里
lr_file <- file.path(MODEL_DIR, "lr_network_human_21122021.rds")     # 不变
if (!file.exists(lt_file)) download.file(
  "https://zenodo.org/record/3260758/files/ligand_target_matrix_nsga2r_final.rds", lt_file, mode = "wb")
if (!file.exists(lr_file)) download.file(
  "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds", lr_file, mode = "wb")
ligand_target_matrix <- readRDS(lt_file)
lr_network <- readRDS(lr_file)

## ---- 2. 读全细胞对象，定义发送/接收细胞 ----
obj <- readRDS(file.path(BASE_DIR, "results", "06_GSE190452_annotated_full.rds"))
obj <- JoinLayers(obj)
meta <- obj@meta.data
ct_col <- colnames(meta)[grepl("cell_type|celltype", colnames(meta), ignore.case = TRUE)][1]
cat("细胞类型列:", ct_col, "\n"); print(table(meta[[ct_col]]))

## 接收细胞 = 神经元；发送细胞 = 微胶质
recv_label <- grep("neuron|Neuron|NEU", unique(as.character(meta[[ct_col]])), value = TRUE)
send_label <- grep("micro|MG", unique(as.character(meta[[ct_col]])), ignore.case = TRUE, value = TRUE)
cat("接收(神经元):", recv_label, "\n发送(微胶质):", send_label, "\n")

## 分组列（TLE vs Normal）
grp_col <- "group"
cat("分组列:", grp_col, "\n"); print(table(meta[[grp_col]]))
if (length(grp_col) == 0) {  ## 没有现成分组列就从样本名推断
  obj$grp <- ifelse(grepl("TLE", meta$orig.ident), "TLE", "Normal")
  grp_col <- "grp"
} else grp_col <- grp_col[1]

## ---- 3. 神经元中 TLE vs Normal 差异基因 = geneset ----
Idents(obj) <- meta[[ct_col]]

neu <- subset(obj, idents = recv_label)
neu <- JoinLayers(neu)
Idents(neu) <- neu@meta.data[[grp_col]]
deg <- FindMarkers(neu, ident.1 = "TLE", ident.2 = "Control",
                   min.pct = 0.05, logfc.threshold = 0.1, verbose = FALSE)
geneset <- rownames(deg)[deg$p_val_adj < 0.05]
cat("神经元 TLE vs Normal DEG 数:", length(geneset), "\n")
cat("其中铁死亡基因:", intersect(FERRO, geneset), "\n")
MODEL_DIR <- file.path(BASE_DIR, "nichenet_models")
ligand_target_matrix <- readRDS(file.path(MODEL_DIR, "ligand_target_matrix.rds"))
lr_network <- readRDS(file.path(MODEL_DIR, "lr_network_human_21122021.rds"))
## ---- 4. NicheNet 配体活性分析 ----
expressed_recv <- rownames(neu)[rowSums(GetAssayData(neu, layer = "data") > 0) / ncol(neu) > 0.05]
mg <- subset(obj, idents = send_label)
expressed_send <- rownames(mg)[rowSums(GetAssayData(mg, layer = "data") > 0) / ncol(mg) > 0.05]

geneset_use <- intersect(geneset, rownames(ligand_target_matrix))
background  <- intersect(expressed_recv, rownames(ligand_target_matrix))
ligands <- lr_network$from %>% unique() %>%
  intersect(expressed_send) %>%
  intersect(colnames(ligand_target_matrix))
receptors   <- lr_network$to %>% unique() %>% intersect(expressed_recv)

act <- predict_ligand_activities(geneset = geneset_use, background_expressed_genes = background,
                                 ligand_target_matrix = ligand_target_matrix, potential_ligands = ligands)
act <- act %>% arrange(-aupr_corrected)
cat("\n配体活性排名 Top10:\n"); print(head(act, 10))
lg3_rank <- which(act$test_ligand == "LGALS3")
cat("\nLGALS3 排名:", lg3_rank, "/", nrow(act),
    " aupr_corrected =", act$aupr_corrected[lg3_rank],
    " P =", act$pearson[lg3_rank] %>% signif(3), "\n")
write.csv(act, file.path(OUT_DIR, "21_nichenet_ligand_activity.csv"), row.names = FALSE)

## ---- 5. LGALS3 的预测靶基因 ----
top_lig <- act$test_ligand[1:min(20, nrow(act))]
lt_use <- ligand_target_matrix[background, top_lig]
lg3_targets <- ligand_target_matrix[background, "LGALS3"]
lg3_top <- sort(lg3_top <- lg3_targets, decreasing = TRUE)[1:30]
cat("\nLGALS3 Top30 预测靶基因:\n"); print(round(lg3_top, 4))

tg <- data.frame(gene = names(lg3_top), weight = as.numeric(lg3_top),
                 in_geneset = names(lg3_top) %in% geneset_use,
                 ferroptosis = names(lg3_top) %in% FERRO)
write.csv(tg, file.path(OUT_DIR, "21_nichenet_LGALS3_targets.csv"), row.names = FALSE)
cat("Top30 中属于神经元DEG:", sum(tg$in_geneset), " 铁死亡基因:", sum(tg$ferroptosis), "\n")
cat("命中铁死亡:", tg$gene[tg$ferroptosis], "\n")

## 富集检验：LGALS3 Top100 靶基因 × 铁死亡基因
top100 <- names(sort(lg3_targets, decreasing = TRUE))[1:100]
ov <- intersect(top100, FERRO)
enrich <- phyper(length(ov) - 1, length(intersect(FERRO, background)),
                 length(background) - length(intersect(FERRO, background)),
                 100, lower.tail = FALSE)
cat(sprintf("\nLGALS3 Top100 靶基因中铁死亡基因: %s (超几何 P = %.3g)\n",
            paste(ov, collapse = ", "), enrich))

## ---- 6. Figure S8 ----
act_top <- head(act, 15) %>%
  mutate(ligand = factor(test_ligand, levels = rev(test_ligand)),
         hl = test_ligand == "LGALS3")
pA <- ggplot(act_top, aes(aupr_corrected, ligand, fill = hl)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("FALSE" = GRAY, "TRUE" = PURPLE), guide = "none") +
  labs(title = "A  Ligand activity (microglia → neurons)",
       x = "AUPR (corrected)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(text = element_text(color = INK), plot.title = element_text(face = "bold", size = 12),
        plot.background = element_rect(fill = PAPER, color = NA),
        panel.background = element_rect(fill = PAPER, color = NA))

tg_top <- head(tg, 20) %>%
  mutate(gene = factor(gene, levels = rev(gene)))
pB <- ggplot(tg_top, aes(weight, gene, fill = ferroptosis)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("FALSE" = GRAY, "TRUE" = PINK), name = "ferroptosis") +
  labs(title = "B  LGALS3 predicted target genes", x = "Regulatory potential", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(text = element_text(color = INK), plot.title = element_text(face = "bold", size = 12),
        plot.background = element_rect(fill = PAPER, color = NA),
        panel.background = element_rect(fill = PAPER, color = NA),
        legend.position = "top")

library(patchwork)
figS8 <- pA | pB
ggsave(file.path(OUT_DIR, "FigureS8_NicheNet.png"), figS8, width = 11, height = 5.5, dpi = 300, bg = PAPER)
ggsave(file.path(OUT_DIR, "FigureS8_NicheNet.pdf"), figS8, width = 11, height = 5.5, bg = PAPER)
cat("\n已保存 21_nichenet_*.csv 与 FigureS8\n")

## ---- 补丁：假设驱动，强制纳入 LGALS3 ----
ligands_fix <- union(ligands, "LGALS3") %>% intersect(colnames(ligand_target_matrix))
cat("LGALS3 在矩阵中:", "LGALS3" %in% ligands_fix, "\n")

act2 <- predict_ligand_activities(geneset = geneset_use,
                                  background_expressed_genes = background,
                                  ligand_target_matrix = ligand_target_matrix,
                                  potential_ligands = ligands_fix)
act2 <- act2 %>% arrange(-aupr_corrected)
i <- which(act2$test_ligand == "LGALS3")
cat(sprintf("LGALS3 排名: %d/%d  aupr_corrected = %.4f  pearson = %.4f\n",
            i, nrow(act2), act2$aupr_corrected[i], act2$pearson[i]))
write.csv(act2, file.path(OUT_DIR, "21_nichenet_ligand_activity_v2.csv"), row.names = FALSE)

## 富集检验重算（Top100 靶基因 × 铁死亡）
lg3_targets <- ligand_target_matrix[background, "LGALS3"]
top100 <- names(sort(lg3_targets, decreasing = TRUE))[1:100]
ov <- intersect(top100, FERRO)
enrich <- phyper(length(ov) - 1, length(intersect(FERRO, background)),
                 length(background) - length(intersect(FERRO, background)),
                 100, lower.tail = FALSE)
cat(sprintf("Top100 靶基因中铁死亡基因: %s (P = %.3g)\n", paste(ov, collapse = ", "), enrich))
