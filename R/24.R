# ============================================================================
# Script      : 24.R
# Purpose     : Post-transcriptional gatekeeper tests: LGALS3 protein-mRNA discordance, mTORC1/translation-initiation, proteasome/autophagy and exosome-loading pathways.
# Manuscript  : Supplementary analyses (post-transcriptional regulation)
# Input       : GSE168375 expression (data/GSE168375_expr.csv or rebuilt from RAW.tar)
# Output      : results/24_* tables/figures
# Dependencies: dplyr, ggplot2, readr, stringr, tibble, tidyr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## =============================================================================
## 24_posttranscriptional_gatekeeper.R  （自包含完整版 v2）
## 方案第三部分（转录后调控版）：
##   3.1 闸门：LGALS3 "蛋白升、mRNA 不升" 不一致性验证
##   3.3 翻译效率通路：mTORC1 / 翻译起始
##   3.4 蛋白降解通路：蛋白酶体 / 自噬-溶酶体 / 分子伴侣（含 HSPB1）
##   3.5 外泌体装载释放：ESCRT / 溶酶体胞吐 / RAB27
##
## 数据：
##   GSE168375 人皮层【纯化星形胶质细胞】RNA-seq，31 epilepsy vs 12 control
##   （注意：非小胶质细胞、非海马，手稿中须如实定位；分组来自 GEO 官方元数据，
##    内嵌于本脚本，peritumor 6 例剔除）
##   表达矩阵：优先读 data/GSE168375_expr.csv；不存在则自动从
##   data/GSE168375_RAW.tar 构建（需 org.Hs.eg.db）
##
## 依赖：tidyverse + org.Hs.eg.db（仅构建矩阵时需要）
## 输出：results/24_posttranscriptional/
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tibble)
  library(stringr); library(tidyr)
})

## --------------------------- 0. 配置区 --------------------------------------
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "results", "24_posttranscriptional")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(sprintf(...), "\n")

GSE168375_FILE <- file.path(DATA_DIR, "GSE168375_expr.csv")
TAR <- file.path(DATA_DIR, "GSE168375_RAW.tar")

## 蛋白侧已知结果（用于不一致性图）
PROTEIN_LOG2FC <- c(rat_proteome = 2.354)   # Lgals3 大鼠 KA 5.1x

## --------------------------- 1. 通路基因集（人工策展） -----------------------
GENE_SETS <- list(
  ## 3.3 翻译效率
  MTORC1_hallmark_core = c("MTOR","RPTOR","MLST8","RPS6KB1","RPS6KB2","EIF4EBP1","EIF4EBP2",
                           "EIF4E","EIF4G1","EIF4G3","RHEB","AKT1","AKT2","TSC1","TSC2",
                           "RPS6","EIF4A1","EIF4A2","DEPTOR"),
  Translation_initiation = c("EIF1","EIF2S1","EIF2S2","EIF2S3","EIF3A","EIF3B","EIF3C","EIF3D",
                             "EIF3E","EIF4E","EIF4G1","EIF4A1","EIF4B","EIF5","PABPC1","PAIP1"),
  ## 3.4 蛋白降解
  Proteasome = c("PSMA1","PSMA2","PSMA3","PSMA4","PSMB1","PSMB2","PSMB3","PSMB4","PSMB5",
                 "PSMC1","PSMC2","PSMD1","PSMD2","PSMD4","PSMD11","PSME1","PSME2","SEM1"),
  Autophagy_lysosome = c("ATG5","ATG7","ATG12","ATG16L1","MAP1LC3B","GABARAP","SQSTM1","BECN1",
                         "ULK1","WIPI2","LAMP1","LAMP2","CTSD","CTSB","ATP6V1A","ATP6V0A1","MCOLN1"),
  Chaperone = c("HSPB1","HSPA1A","HSPA5","HSPA8","HSPA9","HSP90AA1","HSP90AB1","HSPD1",
                "DNAJA1","DNAJA2","DNAJB1","CCT2","CCT3","CCT5","TCP1","PPIA","BAG3","STUB1"),
  ## 3.5 外泌体装载与释放
  ESCRT = c("TSG101","HGS","STAM1","STAM2","CHMP2A","CHMP2B","CHMP4A","CHMP4B","CHMP4C",
            "VPS4A","VPS4B","VTA1","PDCD6IP","SDCBP","MVB12A","IST1"),
  Lysosomal_exocytosis = c("VAMP7","VAMP8","SNAP23","STX4","SYT7","MCOLN1","ARL8A","ARL8B",
                           "BORCS5","RAB3A","RAB27A","RAB27B","LAMP1","ATP6V1H"),
  Rab27_exosome_release = c("RAB27A","RAB27B","SYTL1","SYTL2","SYTL4","EXPH5","MYRIP","UNC13D"),
  ## 对照：核糖体（排除"整体转录活跃"混淆——若它也显著，其他模块的显著须按全局激活解读）
  CTRL_ribosome = c("RPL3","RPL5","RPL7","RPL10","RPS3","RPS6","RPS9","RPS14","RPLP0","RPLP1")
)

## --------------------------- 2. GSM -> 诊断（GEO 官方元数据） -----------------
gsm_map <- tibble::tribble(
  ~gsm, ~dx,
  "GSM5137502","peritumor","GSM5137503","peritumor","GSM5137504","epilepsy",
  "GSM5137505","control","GSM5137506","epilepsy","GSM5137507","epilepsy",
  "GSM5137508","epilepsy","GSM5137509","epilepsy","GSM5137510","epilepsy",
  "GSM5137511","control","GSM5137512","epilepsy","GSM5137513","peritumor",
  "GSM5137514","control","GSM5137515","control","GSM5137516","control",
  "GSM5137517","epilepsy","GSM5137518","epilepsy","GSM5137519","control",
  "GSM5137520","control","GSM5137521","control","GSM5137522","epilepsy",
  "GSM5137523","epilepsy","GSM5137524","epilepsy","GSM5137525","epilepsy",
  "GSM5137526","control","GSM5137527","epilepsy","GSM5137528","epilepsy",
  "GSM5137529","epilepsy","GSM5137530","epilepsy","GSM5137531","peritumor",
  "GSM5137532","control","GSM5137533","epilepsy","GSM5137534","epilepsy",
  "GSM5137535","epilepsy","GSM5137536","peritumor","GSM5137537","control",
  "GSM5137538","epilepsy","GSM5137539","epilepsy","GSM5137540","epilepsy",
  "GSM5137541","epilepsy","GSM5137542","epilepsy","GSM5137543","epilepsy",
  "GSM5137544","epilepsy","GSM5137545","control","GSM5137546","peritumor",
  "GSM5137547","epilepsy","GSM5137548","epilepsy","GSM5137549","epilepsy",
  "GSM5137550","epilepsy"
)

## --------------------------- 3. 构建/读取表达矩阵 ----------------------------
build_from_tar <- function() {
  msg("从 RAW.tar 构建矩阵（仅需一次，之后自动读 csv）……")
  EXD <- file.path(DATA_DIR, "GSE168375_RAW")
  if (!dir.exists(EXD)) untar(TAR, exdir = EXD)
  fs <- list.files(EXD, pattern = "counts\\.txt\\.gz$", full.names = TRUE)
  cts <- lapply(fs, function(f) {
    d <- readr::read_tsv(f, col_names = c("ensg", "count"), show_col_types = FALSE)
    d %>% dplyr::mutate(gsm = str_extract(basename(f), "GSM\\d+"))
  }) %>% bind_rows() %>% dplyr::left_join(gsm_map, by = "gsm")
  cts <- cts %>% dplyr::filter(dx %in% c("epilepsy", "control"))
  mat_ens <- cts %>%
    dplyr::group_by(ensg, gsm) %>%
    dplyr::summarise(count = sum(count), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = gsm, values_from = count, values_fill = 0)
  m <- as.matrix(mat_ens[, -1]); rownames(m) <- mat_ens$ensg
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("需要 org.Hs.eg.db：BiocManager::install('org.Hs.eg.db')")
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(m),
                               column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  keep <- !is.na(sym)
  m2 <- m[keep, ]; rownames(m2) <- toupper(sym[keep])
  m2 <- rowsum(m2, group = rownames(m2))
  write_csv(as_tibble(m2, rownames = "gene"), GSE168375_FILE)
  msg("已存 %s", GSE168375_FILE)
  m2
}

if (file.exists(GSE168375_FILE)) {
  df <- read_csv(GSE168375_FILE, show_col_types = FALSE)
  m2 <- as.matrix(df[, -1]); rownames(m2) <- df$gene
} else {
  if (!file.exists(TAR)) stop("既无 ", GSE168375_FILE, " 也无 ", TAR)
  m2 <- build_from_tar()
}
msg("GSE168375 矩阵：%d 基因 x %d 样本", nrow(m2), ncol(m2))

grp <- ifelse(gsm_map$dx[match(colnames(m2), gsm_map$gsm)] == "epilepsy", "TLE", "Con")
names(grp) <- colnames(m2)
msg("分组：TLE=%d Con=%d（应为 31 vs 12）", sum(grp=="TLE"), sum(grp=="Con"))
stopifnot(sum(grp=="TLE") == 31, sum(grp=="Con") == 12)

## --------------------------- 4. 分析函数 ------------------------------------
one_gene_de <- function(mat, grp, gene = "LGALS3", label = "") {
  if (!gene %in% rownames(mat)) { msg("[%s] !! 矩阵中无 %s", label, gene); return(NULL) }
  x <- as.numeric(mat[gene, ]); l2 <- log2(x + 1)
  fc <- mean(l2[grp == "TLE"]) - mean(l2[grp == "Con"])
  p  <- tryCatch(wilcox.test(l2[grp == "TLE"], l2[grp == "Con"])$p.value, error = function(e) NA)
  tibble(dataset = label, gene = gene, log2FC_mRNA = round(fc, 3),
         wilcox_p = signif(p, 3), detect_rate = round(mean(x > 0), 3),
         n_TLE = sum(grp == "TLE"), n_Con = sum(grp == "Con"))
}

module_scores <- function(mat, grp, label = "") {
  z <- t(scale(t(log2(mat + 1))))
  lapply(names(GENE_SETS), function(gs) {
    g <- intersect(GENE_SETS[[gs]], rownames(z))
    if (length(g) < 5) return(NULL)
    sc <- colMeans(z[g, , drop = FALSE])
    p <- tryCatch(wilcox.test(sc[grp == "TLE"], sc[grp == "Con"])$p.value, error = function(e) NA)
    d <- mean(sc[grp == "TLE"]) - mean(sc[grp == "Con"])
    cor_lg <- if ("LGALS3" %in% rownames(z))
      suppressWarnings(cor(sc, z["LGALS3", ], method = "spearman")) else NA_real_
    tibble(dataset = label, pathway = gs, n_genes = length(g),
           delta_TLE_vs_Con = round(d, 3), wilcox_p = signif(p, 3),
           cor_LGALS3 = round(cor_lg, 3))
  }) %>% bind_rows()
}

## --------------------------- 5. 主流程 --------------------------------------
de_tab <- one_gene_de(m2, grp, "LGALS3", "GSE168375")
sc_tab <- module_scores(m2, grp, "GSE168375")

write_csv(de_tab, file.path(OUT_DIR, "24_gatekeeper_LGALS3_mRNA.csv"))
write_csv(sc_tab, file.path(OUT_DIR, "24_pathway_scores.csv"))
print(de_tab); print(as.data.frame(sc_tab))

## ---- 3.1 闸门判定 ----
msg("===== 闸门判定 =====")
r <- de_tab[1, ]
verdict <- if (abs(r$log2FC_mRNA) < 0.5) {
  "支持转录后主导（mRNA 基本不变）"
} else if (r$log2FC_mRNA < PROTEIN_LOG2FC[1] - 1) {
  "mRNA 升幅远小于蛋白（5.1x），仍以转录后为主"
} else {
  "!! mRNA 显著上调，需改写为 转录+转录后共同作用"
}
msg("  GSE168375: mRNA log2FC=%.2f（P=%.3g；蛋白 log2FC=2.35）→ %s",
    r$log2FC_mRNA, r$wilcox_p, verdict)

## ---- 全局激活检查（核糖体对照） ----
ribo_p <- sc_tab$wilcox_p[sc_tab$pathway == "CTRL_ribosome"]
if (length(ribo_p) > 0 && ribo_p < 0.05) {
  msg("!! 注意：核糖体对照模块也显著（P=%.3g）——癫痫组存在整体转录/翻译活跃，", ribo_p)
  msg("!! 其他模块的显著性应按全局激活解读，不能作为 LGALS3 特异性转录后证据。")
}

## ---- 图 1：mRNA vs 蛋白不一致性 ----
p1 <- ggplot(de_tab, aes(x = dataset, y = log2FC_mRNA)) +
  geom_col(fill = "#0E3A6B", width = 0.55) +
  geom_hline(yintercept = PROTEIN_LOG2FC[1], color = "#C9828B", linetype = 2, linewidth = 1) +
  annotate("text", x = 0.6, y = PROTEIN_LOG2FC[1], label = "Protein: 5.1x (rat KA)",
           hjust = 0, vjust = -0.6, color = "#C9828B", size = 3.6) +
  geom_hline(yintercept = 0, color = "grey60") +
  labs(title = "LGALS3: mRNA vs protein discordance (gatekeeper 3.1)",
       y = "log2FC (TLE vs Con)", x = NULL) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "24_gatekeeper_discordance.pdf"), p1, width = 5.5, height = 4.2)

## ---- 图 2：通路评分 ----
sc_tab$pathway <- factor(sc_tab$pathway, levels = names(GENE_SETS))
p2 <- ggplot(sc_tab, aes(x = pathway, y = delta_TLE_vs_Con, fill = wilcox_p < 0.05)) +
  geom_col(width = 0.65) + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#0E3A6B", "FALSE" = "#B9C6D6"), name = "P < 0.05") +
  geom_hline(yintercept = 0, color = "grey60") +
  labs(title = "Post-transcriptional pathway activity (TLE vs Con, GSE168375)",
       subtitle = "Note: CTRL_ribosome significant => global activation, interpret others with caution",
       y = "Delta module score", x = NULL) +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "24_pathway_scores.pdf"), p2, width = 8, height = 5)

msg("===== 完成，输出目录：%s =====", OUT_DIR)

## ---- GSE88992 修正版：数据已是 log2 尺度，直接算均值差 ----
de92 <- lapply(c("6h","12h","24h"), function(tp) {
  ka <- grep(paste0("kainate_", tp), colnames(m92))
  sa <- grep(paste0("saline_", tp),  colnames(m92))
  x <- m92["LGALS3", ]
  fc <- mean(x[ka]) - mean(x[sa])
  p  <- tryCatch(wilcox.test(x[ka], x[sa])$p.value, error = function(e) NA)
  tibble(dataset = paste0("GSE88992_", tp), gene = "LGALS3",
         log2FC_mRNA = round(fc, 3), wilcox_p = signif(p, 3),
         detect_rate = 1, n_TLE = length(ka), n_Con = length(sa))
}) %>% bind_rows()
print(de92)

## 更新总表与图（删掉旧的错误行再加新行）
de_tab2 <- de_tab2 %>% dplyr::filter(!str_detect(dataset, "GSE88992")) %>% bind_rows(de92)
write_csv(de_tab2, file.path(OUT_DIR, "24_gatekeeper_LGALS3_mRNA.csv"))

p1b <- ggplot(de_tab2, aes(x = dataset, y = log2FC_mRNA)) +
  geom_col(fill = c("#0E3A6B", rep("#7A9CC4", 3)), width = 0.6) +
  geom_hline(yintercept = PROTEIN_LOG2FC[1], color = "#C9828B", linetype = 2, linewidth = 1) +
  annotate("text", x = 0.55, y = PROTEIN_LOG2FC[1], label = "Protein: 5.1x (rat KA, chronic)",
           hjust = 0, vjust = -0.6, color = "#C9828B", size = 3.4) +
  geom_hline(yintercept = 0, color = "grey60") +
  labs(title = "LGALS3 mRNA: chronic human (stable) vs acute mouse (induced)",
       subtitle = "Gatekeeper 3.1: post-transcriptional maintenance is a chronic-phase phenomenon",
       y = "log2FC (TLE/KA vs Con)", x = NULL) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "24_gatekeeper_discordance.pdf"), p1b, width = 6.5, height = 4.2)
msg("已更新（修正 log 尺度）")
