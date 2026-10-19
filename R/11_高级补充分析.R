# ============================================================================
# Script      : 11_高级补充分析.R
# Purpose     : Advanced supplementary single-cell follow-ups: pseudotime checks, regulons (decoupleR), consensus clustering, miRNA links.
# Manuscript  : Supplementary analyses
# Input       : results/06_GSE190452_annotated_full.rds and related intermediates
# Output      : supplementary tables/figures (exploratory)
# Dependencies: CellChat, ConsensusClusterPlus, GSVA, Seurat, SeuratWrappers, data.table, decoupleR, dplyr, ggplot2, ggraph, igraph, monocle3, multiMiR, pROC, tidyr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================

# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# BiocManager::install(c("monocle3", "decoupleR", "ConsensusClusterPlus"),
#                      site_repository = "https://mirrors.westlake.edu.cn/bioconductor")
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# remotes::install_github("satijalab/seurat-wrappers")   # 拟时序需要
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# install.packages(c("ggraph", "igraph"))
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root

obj <- readRDS(file.path(BASE_DIR, "results/06_GSE190452_annotated_full.rds"))
grep("LDAM|Ferro|UCell|score", colnames(obj@meta.data), value = TRUE, ignore.case = TRUE)
dfp <- mg@meta.data[!is.na(mg$pseudotime), ]
# LDAM 评分是否沿 pseudotime 升高？（验证根方向）
cor.test(dfp$pseudotime, dfp$LDAM_UCell, method = "spearman")
# LDAM 与铁死亡直接相关
cor.test(dfp$LDAM_UCell, dfp$Ferroptosis_UCell, method = "spearman")
# 分组看：LDAM-high vs LDAM-low 的铁死亡差异
wilcox.test(Ferroptosis_UCell ~ LDAM_level, data = dfp)
tapply(dfp$Ferroptosis_UCell, dfp$LDAM_level, mean, na.rm = TRUE)
## ============================================================
## 11_高级补充分析.R
## 补1: 拟时序(monocle3)      补2: 转录因子(decoupleR/DoRothEA)
## 补3: LOOCV + 逐基因ROC      补4: CMap药物预测签名文件
## 补5: 铁死亡分型(Consensus)  补6: 血液LDAM样单核评分
## 补7: 免疫浸润(修复版)       补8: miRNA调控网络
## 补9: 各细胞类型LDAM评分热图
## ============================================================
library(CellChat)
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR <- file.path(BASE_DIR, "results"); FIG_DIR <- file.path(BASE_DIR, "figures")
ADV_DIR <- file.path(FIG_DIR, "advanced"); dir.create(ADV_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr)
  library(Seurat); library(pROC)
})
COL <- c(Control = "#7BA7BC", TLE = "#C0392B", Epilepsy = "#C0392B",
         `LDAM-low` = "#7BA7BC", `LDAM-high` = "#C0392B")
hub_genes <- readLines(file.path(OUT_DIR, "04_hub_genes.txt"))
logit <- function(...) cat(..., "\n")

## ============================================================
## 补 9｜各细胞类型 LDAM / 铁死亡评分（先跑这个，补1/2都依赖 obj）
## ============================================================
obj <- mg <- NULL
tryCatch({
  obj <- readRDS(file.path(OUT_DIR, "06_GSE190452_annotated_full.rds"))
  mg  <- readRDS(file.path(OUT_DIR, "06_GSE190452_microglia.rds"))
  mg <- JoinLayers(mg)
  ct_scores <- obj@meta.data %>%
    group_by(cell_type) %>%
    summarise(LDAM = mean(LDAM_UCell, na.rm = TRUE),
              Ferroptosis = mean(Ferroptosis_UCell, na.rm = TRUE))
  ct_l <- pivot_longer(ct_scores, c(LDAM, Ferroptosis), names_to = "pathway", values_to = "score")
  p9 <- ggplot(ct_l, aes(reorder(cell_type, score), score, fill = pathway)) +
    geom_col(position = "dodge", width = 0.65) +
    scale_fill_manual(values = c(LDAM = "#C0392B", Ferroptosis = "#5B8E7D")) +
    coord_flip() + theme_classic(base_size = 10) +
    labs(title = "LDAM & Ferroptosis scores across cell types", x = NULL)
  ggsave(file.path(ADV_DIR, "11_9_celltype_pathway_scores.pdf"), p9, width = 6.5, height = 4)
  fwrite(ct_scores, file.path(OUT_DIR, "11_9_celltype_pathway_scores.csv"))
  logit("补9 完成")
}, error = function(e) logit("补9 失败: ", e$message))

## ============================================================
## 补 1｜拟时序分析（monocle3，LDAM-low → LDAM-high 轨迹）
## ============================================================
tryCatch({
  if (is.null(mg)) stop("mg 未加载")
  library(monocle3); library(SeuratWrappers)
  cds <- as.cell_data_set(mg)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  ## 以 LDAM-low 组为起点
  root_cells <- colnames(cds)[mg$LDAM_level == "LDAM-low" & !is.na(mg$LDAM_level)]
  cds <- order_cells(cds, root_cells = root_cells)
  pt <- pseudotime(cds)
  mg$pseudotime <- pt[colnames(mg)]
  pdf(file.path(ADV_DIR, "11_1_pseudotime_umap.pdf"), width = 12, height = 5)
  print(plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE) + ggtitle("Pseudotime"))
  print(plot_cells(cds, color_cells_by = "LDAM_level", label_cell_groups = FALSE,
                   label_leaves = FALSE, label_branch_points = FALSE) + ggtitle("LDAM level"))
  dev.off()
  ## 关键结果：铁死亡评分沿轨迹升高
  dfp <- mg@meta.data[!is.na(mg$pseudotime), ]
  ct <- cor.test(dfp$pseudotime, dfp$Ferroptosis_UCell, method = "spearman")
  p1c <- ggplot(dfp, aes(pseudotime, Ferroptosis_UCell)) +
    geom_point(alpha = 0.2, size = 0.4, color = "#5B8E7D") +
    geom_smooth(method = "loess", color = "#C0392B") +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
             label = paste0("rho = ", round(ct$estimate, 3), ", P ", format(ct$p.value, digits = 2))) +
    theme_classic(base_size = 10) +
    labs(title = "Ferroptosis score along LDAM trajectory")
  ggsave(file.path(ADV_DIR, "11_1_ferroptosis_along_pseudotime.pdf"), p1c, width = 5.5, height = 4.5)
  fwrite(data.frame(stat = c("pseudotime_vs_ferroptosis_rho", "p"),
                    value = c(ct$estimate, ct$p.value)),
         file.path(OUT_DIR, "11_1_pseudotime_correlation.csv"))
  logit("补1 完成: ferroptosis ~ pseudotime rho =", round(ct$estimate, 3),
        ", P =", format(ct$p.value, digits = 2))
}, error = function(e) logit("补1 失败(多为monocle3/SeuratWrappers未装): ", e$message))

## ============================================================
## 补 2｜转录因子活性分析（decoupleR + DoRothEA）
## ============================================================
tryCatch({
  if (is.null(mg)) stop("mg 未加载")
  library(decoupleR)
  ## LDAM-high vs low 的 DE t 值作为输入
  Idents(mg) <- "LDAM_level"
  de <- FindMarkers(mg, ident.1 = "LDAM-high", ident.2 = "LDAM-low",
                    min.pct = 0.05, logfc.threshold = 0)
  tvec <- de$avg_log2FC * -log10(de$p_val + 1e-300)   # 近似 t 统计量
  names(tvec) <- rownames(de)
  net <- decoupleR::get_dorothea(organism = "human", levels = c("A", "B", "C"))
  tf_res <- run_ulm(mat = matrix(tvec, ncol = 1, dimnames = list(names(tvec), "LDAMhigh_vs_low")),
                    net = net, .source = "tf", .target = "target", minsize = 5)
  fwrite(tf_res, file.path(OUT_DIR, "11_2_TF_activity_LDAMhigh.csv"))
  top_tf <- tf_res %>% arrange(p_value) %>% head(20)
  top_tf$tf <- factor(top_tf$tf, levels = rev(top_tf$tf))
  p2 <- ggplot(top_tf, aes(tf, score, fill = score > 0)) +
    geom_col(width = 0.6) + coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2C5F8A"), guide = "none") +
    theme_classic(base_size = 10) +
    labs(title = "Top TF activities: LDAM-high vs LDAM-low", x = NULL, y = "ULM score")
  ggsave(file.path(ADV_DIR, "11_2_TF_activity.pdf"), p2, width = 6, height = 6)
  logit("补2 完成, top TF: ", paste(head(top_tf$tf, 8), collapse = ", "))
}, error = function(e) logit("补2 失败: ", e$message))

## ============================================================
## 补 3｜LOOCV + 逐基因 ROC + Bootstrap CI（发现队列 GSE168375）
## ============================================================
tryCatch({
  ex <- fread(file.path(OUT_DIR, "02_GSE168375_vst_normalized.csv")) |>
    tibble::column_to_rownames("gene") |> as.matrix()
  mt <- fread(file.path(OUT_DIR, "02_GSE168375_sample_meta.csv"))
  sm <- intersect(mt$sample, colnames(ex))
  grp <- factor(ifelse(mt$group[match(sm, mt$sample)] %in% c("TLE","Epilepsy"), "Epilepsy", "Control"),
                levels = c("Control", "Epilepsy"))
  X <- t(ex[hub_genes, sm])
  ## 3.1 LOOCV
  pred <- numeric(length(grp))
  for (i in seq_along(grp)) {
    fit <- glm(grp[-i] ~ ., data = data.frame(grp = grp[-i], X[-i, ]), family = binomial)
    pred[i] <- predict(fit, newdata = data.frame(X[i, , drop = FALSE]), type = "response")
  }
  roc_loocv <- roc(grp, pred, quiet = TRUE); ci_l <- ci.auc(roc_loocv)
  logit("补3 LOOCV AUC =", round(auc(roc_loocv), 3),
        " (95%CI ", round(ci_l[1], 2), "-", round(ci_l[3], 2), ")")
  ## 3.2 逐基因 ROC
  roc_list <- lapply(hub_genes, function(g) roc(grp, X[, g], quiet = TRUE))
  names(roc_list) <- hub_genes
  pdf(file.path(ADV_DIR, "11_3_LOOCV_pergene_ROC.pdf"), width = 11, height = 5)
  par(mfrow = c(1, 2))
  plot(roc_loocv, col = "#2C5F8A", lwd = 2,
       main = paste0("LOOCV ROC, AUC = ", round(auc(roc_loocv), 2)))
  abline(0, 1, lty = 2, col = "grey60")
  cols4 <- c("#C0392B", "#5B8E7D", "#D4A24E", "#2C5F8A")
  plot(roc_list[[1]], col = cols4[1], lwd = 2, main = "Per-gene ROC")
  for (k in 2:4) plot(roc_list[[k]], col = cols4[k], lwd = 2, add = TRUE)
  legend("bottomright", paste(hub_genes, round(sapply(roc_list, auc), 2)),
         col = cols4, lwd = 2, cex = 0.8)
  abline(0, 1, lty = 2, col = "grey60"); dev.off()
  ## 3.3 Bootstrap CI（1000 次重抽样）
  set.seed(42)
  aucs <- replicate(1000, {
    idx <- sample(length(grp), replace = TRUE)
    if (length(unique(grp[idx])) < 2) return(NA)
    fit <- tryCatch(glm(grp[idx] ~ ., data = data.frame(grp = grp[idx], X[idx, ]), family = binomial),
                    error = function(e) NULL)
    if (is.null(fit)) return(NA)
    auc(roc(grp[idx], predict(fit, type = "response"), quiet = TRUE))
  })
  aucs <- aucs[!is.na(aucs)]
  logit("补3 Bootstrap AUC =", round(mean(aucs), 3),
        " (95%CI ", round(quantile(aucs, 0.025), 2), "-", round(quantile(aucs, 0.975), 2), ")")
  fwrite(data.frame(metric = c("LOOCV_AUC","LOOCV_CI_lo","LOOCV_CI_hi","Boot_AUC","Boot_CI_lo","Boot_CI_hi"),
                    value = c(auc(roc_loocv), ci_l[1], ci_l[3], mean(aucs),
                              quantile(aucs, 0.025), quantile(aucs, 0.975))),
         file.path(OUT_DIR, "11_3_ROC_robustness.csv"))
  logit("补3 完成")
}, error = function(e) logit("补3 失败: ", e$message))

## ============================================================
## 补 4｜CMap/LINCS 药物预测：生成签名文件（上传 clue.io 查询）
## ============================================================
tryCatch({
  fde <- list.files(OUT_DIR, pattern = "09.*LDAM.*(DE|marker).*csv", full.names = TRUE)
  if (length(fde) == 0) fde <- list.files(OUT_DIR, pattern = "LDAMhigh", full.names = TRUE)
  de <- fread(fde[1])
  gcol <- intersect(names(de), c("gene", "Gene", "symbol"))[1]
  fccol <- intersect(names(de), c("avg_log2FC", "log2FoldChange"))[1]
  de <- de[order(-de[[fccol]]), ]
  up150   <- na.omit(de[[gcol]][de[[fccol]] > 0])[1:150]
  down150 <- na.omit(de[[gcol]][de[[fccol]] < 0])[1:150]
  writeLines(up150,   file.path(OUT_DIR, "11_4_CMap_signature_UP150.txt"))
  writeLines(down150, file.path(OUT_DIR, "11_4_CMap_signature_DOWN150.txt"))
  logit("补4 完成: 签名文件已生成。上传到 https://clue.io/query 查询可逆转 LDAM-high 表型的小分子")
}, error = function(e) logit("补4 失败: ", e$message))

## ============================================================
## 补 5｜铁死亡亚型分型（ConsensusClusterPlus，GSE63808）
## ============================================================
tryCatch({
  library(ConsensusClusterPlus)
  ex63 <- fread(file.path(OUT_DIR, "03_GSE63808_normalized.csv")) |>
    tibble::column_to_rownames(1) |> as.matrix()
  ferro <- readLines(file.path(BASE_DIR, "genesets", "Ferroptosis_geneset.txt"))
  fg <- intersect(ferro, rownames(ex63))
  m <- ex63[fg, ]; m <- m[apply(m, 1, sd) > 0, ]
  m <- sweep(m, 1, rowMeans(m))  # 中心化
  cc <- ConsensusClusterPlus(m, maxK = 4, reps = 200, pItem = 0.8,
                             clusterAlg = "hc", distance = "spearman",
                             title = file.path(ADV_DIR, "11_5_consensus"))
  sub <- cc[[2]]$consensusClass   # 默认 2 亚型
  ## 亚型间 LDAM 评分比较
  library(GSVA)
  ldam <- readLines(file.path(BASE_DIR, "genesets", "LDAM_geneset.txt"))
  par <- gsvaParam(ex63, list(LDAM = intersect(ldam, rownames(ex63)),
                              Ferroptosis = fg))
  sc63 <- as.data.frame(t(gsva(par)))
  sc63$subtype <- factor(paste0("C", sub[match(rownames(sc63), colnames(ex63))]))
  wt1 <- wilcox.test(LDAM ~ subtype, sc63); wt2 <- wilcox.test(Ferroptosis ~ subtype, sc63)
  sc_l <- pivot_longer(sc63, c(LDAM, Ferroptosis), names_to = "pathway", values_to = "score")
  p5 <- ggplot(sc_l, aes(subtype, score, fill = subtype)) +
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~pathway, scales = "free") +
    scale_fill_manual(values = c(C1 = "#7BA7BC", C2 = "#C0392B")) +
    theme_classic(base_size = 10) +
    labs(title = paste0("Ferroptosis subtypes (GSE63808): LDAM P = ", format(wt1$p.value, digits = 2),
                        ", Ferroptosis P = ", format(wt2$p.value, digits = 2)))
  ggsave(file.path(ADV_DIR, "11_5_subtype_scores.pdf"), p5, width = 7, height = 4.5)
  fwrite(sc63, file.path(OUT_DIR, "11_5_ferroptosis_subtypes.csv"), row.names = TRUE)
  logit("补5 完成: LDAM P =", format(wt1$p.value, digits = 2),
        ", Ferroptosis P =", format(wt2$p.value, digits = 2))
}, error = function(e) logit("补5 失败(检查03_GSE63808_normalized.csv列名): ", e$message))

## ============================================================
## 补 6｜血液 LDAM 样单核细胞评分（LDAM-high 签名 → 血液 GSVA）
## ============================================================
tryCatch({
  fde <- list.files(OUT_DIR, pattern = "09.*LDAM.*(DE|marker).*csv", full.names = TRUE)
  de <- fread(fde[1])
  gcol <- intersect(names(de), c("gene", "Gene", "symbol"))[1]
  fccol <- intersect(names(de), c("avg_log2FC", "log2FoldChange"))[1]
  sig <- de[[gcol]][de[[fccol]] > 0.25][1:100]   # LDAM-high 上调 top100
  eb <- fread(file.path(OUT_DIR, "08_GSE143272_normalized.csv")) |>
    tibble::column_to_rownames(1) |> as.matrix()
  mb <- fread(file.path(OUT_DIR, "08_GSE143272_sample_meta.csv"))
  par_b <- gsvaParam(eb, list(LDAMhigh_signature = intersect(sig, rownames(eb))))
  sb <- data.frame(sample = colnames(eb), score = as.numeric(gsva(par_b)[1, ]))
  sb$group <- mb$group[match(sb$sample, mb$sample)]
  wt <- wilcox.test(score ~ group, sb)
  p6 <- ggplot(sb, aes(group, score, fill = group)) +
    geom_boxplot(outlier.size = 0.5) + geom_point(alpha = 0.4, size = 0.8) +
    scale_fill_manual(values = COL) + theme_classic(base_size = 10) +
    labs(title = paste0("Blood LDAM-high signature score, P = ", format(wt$p.value, digits = 2)),
         y = "GSVA score") + theme(legend.position = "none")
  ggsave(file.path(ADV_DIR, "11_6_blood_LDAMsignature.pdf"), p6, width = 4.5, height = 4.5)
  fwrite(sb, file.path(OUT_DIR, "11_6_blood_LDAMsignature_scores.csv"))
  logit("补6 完成: 血液 LDAM-high 签名 P =", format(wt$p.value, digits = 2))
}, error = function(e) logit("补6 失败: ", e$message))

## ============================================================
## 补 7｜免疫浸润箱线图（修复分组 bug 版）
## ============================================================
tryCatch({
  imm <- fread(file.path(OUT_DIR, "05_ssGSEA_immune_scores.csv"))
  mt1 <- fread(file.path(OUT_DIR, "01_GSE134697_sample_meta.csv"))
  il <- melt(imm, id.vars = "cell_type", variable.name = "sample", value.name = "score")
  il$sample <- gsub("\\.", "-", il$sample)
  il$group <- mt1$tissue[match(il$sample, mt1$sample)]
  il <- il[!is.na(group)]
  pvals <- il[, { w <- wilcox.test(score ~ group)
  .(p = w$p.value, y = max(score) + 0.1*(max(score)-min(score))) }, by = cell_type]
  pvals$label <- ifelse(pvals$p < 0.001, "***", ifelse(pvals$p < 0.01, "**",
                                                       ifelse(pvals$p < 0.05, "*", "ns")))
  p7 <- ggplot(il, aes(group, score, fill = group)) +
    geom_boxplot(outlier.size = 0.4, width = 0.6) +
    facet_wrap(~cell_type, scales = "free", nrow = 2) +
    scale_fill_manual(values = c(Neocortex = "#7BA7BC", Hippocampus = "#C0392B")) +
    geom_text(data = pvals, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 3.5) +
    theme_classic(base_size = 10) +
    labs(title = "Immune infiltration (ssGSEA, GSE134697)", x = NULL) +
    theme(legend.position = "none")
  ggsave(file.path(ADV_DIR, "11_7_immune_infiltration_fixed.pdf"), p7, width = 9, height = 6)
  logit("补7 完成, 显著细胞: ", paste(pvals$cell_type[pvals$p < 0.05], collapse = ", "))
}, error = function(e) logit("补7 失败: ", e$message))

## ============================================================
## 补 8｜hub 基因 miRNA 调控网络
## ============================================================
tryCatch({
  ## 用 multiMiR 查已验证的 miRNA-靶基因互作
  if (!requireNamespace("multiMiR", quietly = TRUE))
    # [setup] one-time package installation commented out for submission;
    # [setup] run 00_setup_genesets.R (or install manually) before first use:
#     BiocManager::install("multiMiR", site_repository = "https://mirrors.westlake.edu.cn/bioconductor", update = FALSE, ask = FALSE)
  library(multiMiR); library(igraph); library(ggraph)
  mm <- get_multimir(org = "hsa", target = hub_genes, table = "validated", summary = TRUE)
  edges <- unique(mm@data[, c("mature_mirna_id", "target_symbol")])
  edges <- edges[edges$mature_mirna_id != "", ]
  fwrite(edges, file.path(OUT_DIR, "11_8_miRNA_hub_network.csv"))
  g <- graph_from_data_frame(edges, directed = FALSE)
  V(g)$type <- V(g)$name %in% hub_genes
  pdf(file.path(ADV_DIR, "11_8_miRNA_network.pdf"), width = 8, height = 8)
  print(ggraph(g, layout = "fr") +
          geom_edge_link(alpha = 0.3, color = "grey60") +
          geom_node_point(aes(color = type, size = type)) +
          scale_color_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2C5F8A"), guide = "none") +
          scale_size_manual(values = c("TRUE" = 5, "FALSE" = 2), guide = "none") +
          geom_node_text(aes(label = ifelse(type, name, "")), repel = TRUE, fontface = "bold") +
          theme_void() + ggtitle("miRNA - hub gene network (validated)"))
  dev.off()
  logit("补8 完成: ", nrow(edges), " 对已验证互作, ", length(unique(edges$mature_mirna_id)), " 个miRNA")
}, error = function(e) logit("补8 失败(网络问题可跳过,不影响其他): ", e$message))

cat("\n== 高级补充分析全部结束，输出目录: ", ADV_DIR, " ==\n")
print(list.files(ADV_DIR))
