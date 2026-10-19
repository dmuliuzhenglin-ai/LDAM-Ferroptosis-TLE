# ============================================================================
# Script      : 糖基化.R
# Purpose     : Validation analysis of the glycosylation-LDAM-ferroptosis axis reusing existing results (offline).
# Manuscript  : Supplementary analyses (glycosylation)
# Input       : existing results/ files
# Output      : figures/glycosylation/ + results tables
# Dependencies: AnnotationDbi, GSVA, Seurat, UCell, data.table, dplyr, ggplot2, org.Hs.eg.db, tidyr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================

## ============================================================
## 12_糖基化验证分析.R
## 验证性分析：糖基化与 LDAM-铁死亡轴的关系
## 数据源：全部复用现有 results/ 文件，不联网
## ============================================================

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR <- file.path(BASE_DIR, "results"); FIG_DIR <- file.path(BASE_DIR, "figures")
GLY_DIR <- file.path(FIG_DIR, "glycosylation"); dir.create(GLY_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr)
})
logit <- function(...) cat(..., "\n")

## ---------- 1. 构建糖基化基因集（org.Hs.eg.db 本地注释，不联网）----------
gly_sets <- NULL
tryCatch({
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    # [setup] one-time package installation commented out for submission;
    # [setup] run 00_setup_genesets.R (or install manually) before first use:
#     BiocManager::install("org.Hs.eg.db",
#                          site_repository = "https://mirrors.westlake.edu.cn/bioconductor",
#                          update = FALSE, ask = FALSE)
  library(org.Hs.eg.db); library(AnnotationDbi)
  go2gene <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = c("GO:0006486",   # protein glycosylation
                                            "GO:0030259",   # lipid glycosylation
                                            "GO:0006487",   # protein N-linked glycosylation
                                            "GO:0001573"),  # ganglioside metabolic process
                                   columns = c("SYMBOL"), keytype = "GO")
  gly_sets <- list(
    Protein_glycosylation = unique(go2gene$SYMBOL[go2gene$GO == "GO:0006486"]),
    Lipid_glycosylation   = unique(go2gene$SYMBOL[go2gene$GO == "GO:0030259"]),
    N_linked_glycosylation= unique(go2gene$SYMBOL[go2gene$GO == "GO:0006487"]),
    Ganglioside_metabolism= unique(go2gene$SYMBOL[go2gene$GO == "GO:0001573"]),
    OGlcNAc_machinery     = c("OGT", "OGT", "OGA", "MGEA5", "GFPT1", "GFPT2",
                              "GNPNAT1", "PGM3", "UAP1", "UAP1L1", "NUP62", "HCFC1"),
    Glycosyltransferases  = unique(go2gene$SYMBOL[grepl("MGAT|B4GALT|B3GALT|ST3GAL|ST6GAL|ST8SIA|FUT", go2gene$SYMBOL)])
  )
  gly_sets <- lapply(gly_sets, na.omit)
  sapply(gly_sets, function(s) logit("基因集: ", s[1], " ... 共 ", length(s), " 个"))
  saveRDS(gly_sets, file.path(OUT_DIR, "12_glycosylation_genesets.rds"))
}, error = function(e) logit("基因集构建失败: ", e$message))

if (is.null(gly_sets)) stop("基因集没建成，把上面的失败信息发我")

## ---------- 2. Bulk 海马（GSE134697）：糖基化评分 + 与 LDAM/铁死亡相关 ----------
tryCatch({
  library(GSVA)
  f1 <- list.files(OUT_DIR, pattern = "01_GSE134697.*(norm|vst|express).*csv", full.names = TRUE)
  logit("bulk 表达矩阵候选: ", paste(basename(f1), collapse = ", "))
  stopifnot(length(f1) > 0)
  ex1 <- fread(f1[1]); gcol <- names(ex1)[1]
  exm <- as.matrix(ex1[, -1, with = FALSE]); rownames(exm) <- ex1[[gcol]]
  sc1 <- as.data.frame(t(gsva(gsvaParam(exm, gly_sets))))
  sc1$sample <- rownames(sc1)
  mt1 <- fread(file.path(OUT_DIR, "01_GSE134697_sample_meta.csv"))
  sc1$group <- mt1$tissue[match(gsub("\\.", "-", sc1$sample), mt1$sample)]
  ## 和已有 LDAM/铁死亡评分合并（05 的 GSVA 结果）
  f5 <- list.files(OUT_DIR, pattern = "05.*GSVA.*score|05_GSVA", full.names = TRUE)
  if (length(f5) > 0) {
    sc5 <- fread(f5[1])
    sc5$sample <- gsub("\\.", "-", sc5[[1]])
    sc1$sample2 <- gsub("\\.", "-", sc1$sample)
    sc1 <- merge(sc1, sc5, by.x = "sample2", by.y = "sample", all.x = TRUE)
  }
  ## 组间比较 + 相关性
  res <- list()
  for (gs in names(gly_sets)) {
    d <- sc1[!is.na(sc1$group), ]
    wt <- wilcox.test(d[[gs]] ~ d$group)
    res[[gs]] <- data.frame(geneset = gs, p = wt$p.value,
                            median_hippo = median(d[[gs]][d$group == "Hippocampus"], na.rm = TRUE),
                            median_ctx   = median(d[[gs]][d$group == "Neocortex"], na.rm = TRUE))
  }
  res_df <- do.call(rbind, res)
  logit("== bulk 组间比较 (Hippocampus vs Neocortex) =="); print(res_df)
  fwrite(res_df, file.path(OUT_DIR, "12_bulk_glycosylation_comparison.csv"))
  ## 与 LDAM/铁死亡评分相关
  if ("LDAM" %in% names(sc1)) {
    cors <- sapply(names(gly_sets), function(gs) {
      ct <- cor.test(sc1[[gs]], sc1$LDAM, method = "spearman")
      c(rho = ct$estimate, p = ct$p.value)
    })
    logit("== 与 LDAM 评分相关性 =="); print(round(t(cors), 4))
  }
  ## 画图
  sc_l <- pivot_longer(sc1[!is.na(sc1$group), ], all_of(names(gly_sets)),
                       names_to = "geneset", values_to = "score")
  p_bulk <- ggplot(sc_l, aes(group, score, fill = group)) +
    geom_boxplot(outlier.size = 0.4) + facet_wrap(~geneset, scales = "free", nrow = 2) +
    scale_fill_manual(values = c(Neocortex = "#7BA7BC", Hippocampus = "#C0392B")) +
    theme_classic(base_size = 9) +
    labs(title = "Glycosylation-related scores: hippocampus vs neocortex") +
    theme(legend.position = "none")
  ggsave(file.path(GLY_DIR, "12_bulk_glycosylation_scores.pdf"), p_bulk, width = 10, height = 6)
  logit("bulk 部分完成")
}, error = function(e) logit("bulk 部分失败(找不到表达矩阵就把候选文件名发我): ", e$message))

## ---------- 3. 单细胞：LDAM-high vs low 糖基化评分 + 相关性 ----------
tryCatch({
  library(Seurat); library(UCell)
  mg <- readRDS(file.path(OUT_DIR, "06_GSE190452_microglia.rds")); mg <- JoinLayers(mg)
  mg <- AddModuleScore_UCell(mg, features = gly_sets, name = "_glyUCell")
  md <- mg@meta.data
  score_cols <- paste0(names(gly_sets), "_glyUCell")
  ## 组间比较
  res2 <- lapply(score_cols, function(sc) {
    wt <- wilcox.test(md[[sc]] ~ md$LDAM_level)
    data.frame(geneset = sc, p = wt$p.value)
  }) |> bind_rows()
  logit("== 单细胞 LDAM-high vs low 糖基化评分 =="); print(res2)
  fwrite(res2, file.path(OUT_DIR, "12_sc_glycosylation_LDAMhigh.csv"))
  ## 与 LDAM/铁死亡评分相关
  cors2 <- sapply(score_cols, function(sc) {
    ct <- cor.test(md[[sc]], md$LDAM_UCell, method = "spearman")
    c(rho = ct$estimate, p = ct$p.value)
  })
  logit("== 单细胞水平与 LDAM 评分相关 =="); print(round(t(cors2), 4))
  ## 关键糖基化酶在 LDAM-high/low 的表达
  key_enzymes <- intersect(c("OGT","OGA","GFPT1","GFPT2","MGAT1","MGAT2","MGAT3","MGAT4A","MGAT5",
                             "B4GALT1","B4GALT5","ST3GAL5","UGCG","HEXA","HEXB","NEU1","GLB1"),
                           rownames(mg))
  Idents(mg) <- "LDAM_level"
  av <- AverageExpression(mg, features = key_enzymes, group.by = "LDAM_level")$RNA
  av_df <- as.data.frame(av) |> tibble::rownames_to_column("gene")
  fwrite(av_df, file.path(OUT_DIR, "12_glyco_enzymes_LDAMlevel.csv"))
  ## 画小提琴图（糖基化评分）
  md_l <- pivot_longer(md, all_of(score_cols), names_to = "geneset", values_to = "score")
  md_l$geneset <- sub("_glyUCell", "", md_l$geneset)
  p_sc <- ggplot(md_l, aes(LDAM_level, score, fill = LDAM_level)) +
    geom_violin(scale = "width") + geom_boxplot(width = 0.1, outlier.size = 0.2) +
    facet_wrap(~geneset, scales = "free", nrow = 2) +
    scale_fill_manual(values = c("LDAM-low" = "#7BA7BC", "LDAM-high" = "#C0392B")) +
    theme_classic(base_size = 9) +
    labs(title = "Glycosylation scores in microglia: LDAM-high vs low") +
    theme(legend.position = "none")
  ggsave(file.path(GLY_DIR, "12_sc_glycosylation_scores.pdf"), p_sc, width = 10, height = 6)
  logit("单细胞部分完成")
}, error = function(e) logit("单细胞部分失败: ", e$message))

## ---------- 4. 外周血（GSE143272）----------
tryCatch({
  library(GSVA)
  eb <- fread(file.path(OUT_DIR, "08_GSE143272_normalized.csv"))
  gcol <- names(eb)[1]
  ebm <- as.matrix(eb[, -1, with = FALSE]); rownames(ebm) <- eb[[gcol]]
  ebm <- ebm[!is.na(rownames(ebm)) & rownames(ebm) != "", ]
  mb <- fread(file.path(OUT_DIR, "08_GSE143272_sample_meta.csv"))
  scb <- as.data.frame(t(gsva(gsvaParam(ebm, gly_sets))))
  scb$sample <- rownames(scb)
  scb$group <- mb$group[match(scb$sample, mb$sample)]
  scb <- scb[!is.na(scb$group), ]
  res3 <- lapply(names(gly_sets), function(gs) {
    wt <- wilcox.test(scb[[gs]] ~ scb$group)
    data.frame(geneset = gs, p = wt$p.value)
  }) |> bind_rows()
  logit("== 血液糖基化评分 (Epilepsy vs Control) =="); print(res3)
  fwrite(res3, file.path(OUT_DIR, "12_blood_glycosylation.csv"))
  scb_l <- pivot_longer(scb, all_of(names(gly_sets)), names_to = "geneset", values_to = "score")
  p_b <- ggplot(scb_l, aes(group, score, fill = group)) +
    geom_boxplot(outlier.size = 0.3) + facet_wrap(~geneset, scales = "free", nrow = 2) +
    scale_fill_manual(values = c(Control = "#7BA7BC", Epilepsy = "#C0392B")) +
    theme_classic(base_size = 9) +
    labs(title = "Glycosylation scores in peripheral blood") +
    theme(legend.position = "none")
  ggsave(file.path(GLY_DIR, "12_blood_glycosylation_scores.pdf"), p_b, width = 10, height = 6)
  logit("血液部分完成")
}, error = function(e) logit("血液部分失败: ", e$message))

## ---------- 5. 小鼠 KA 模型（基因名大小写转换）----------
tryCatch({
  library(GSVA)
  to_mouse <- function(g) paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
  gly_mouse <- lapply(gly_sets, to_mouse)
  f8 <- list.files(OUT_DIR, pattern = "08_GSE88992.*(norm|express).*csv", full.names = TRUE)
  logit("小鼠表达矩阵候选: ", paste(basename(f8), collapse = ", "))
  stopifnot(length(f8) > 0)
  em <- fread(f8[1]); gcol <- names(em)[1]
  emm <- as.matrix(em[, -1, with = FALSE]); rownames(emm) <- em[[gcol]]
  scm <- as.data.frame(t(gsva(gsvaParam(emm, gly_mouse))))
  meta8 <- fread(file.path(OUT_DIR, "08_GSE88992_sample_meta.csv"))
  scm$sample <- rownames(scm)
  scm$group <- meta8$group[match(scm$sample, meta8$sample)]
  res4 <- lapply(names(gly_mouse), function(gs) {
    wt <- tryCatch(wilcox.test(scm[[gs]] ~ scm$group), error = function(e) NULL)
    if (is.null(wt)) return(NULL)
    data.frame(geneset = gs, p = wt$p.value)
  }) |> bind_rows()
  logit("== 小鼠糖基化评分 (KA vs Saline) =="); print(res4)
  fwrite(res4, file.path(OUT_DIR, "12_mouse_glycosylation.csv"))
  logit("小鼠部分完成")
}, error = function(e) logit("小鼠部分失败: ", e$message))

cat("\n== 糖基化验证分析结束, 输出目录: ", GLY_DIR, " ==\n")
print(list.files(GLY_DIR))
