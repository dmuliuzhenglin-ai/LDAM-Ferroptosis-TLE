# ============================================================================
# Script      : 23_exosome_proteome_enrichment.R
# Purpose     : Exosome-protein enrichment in the rat KA hippocampal TMT proteome; Gal-3/CD9 vesicle localization; optional cross-species repeat.
# Manuscript  : Supplementary analyses (cross-species proteome)
# Input       : Differentially_expressed_statistics.xlsx (rat TMT); ExoCarta/Vesiclepedia lists (download or built-in top100)
# Output      : results/23_exosome_proteome/ CSV + PDF
# Dependencies: dplyr, ggplot2, readr, stringr, tibble, tidyr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## =============================================================================
## 23_exosome_proteome_enrichment.R
## 方案4.5：差异蛋白中的外泌体蛋白富集 + Gal-3/CD9 外泌体定位验证 + 跨物种重复
##
## 输入：
##   1) 大鼠 KA 海马 TMT 蛋白质组"全量定量表"（含所有检测到的蛋白 + log2FC + P）
##      —— 注意：背景必须是"所有可检测蛋白"，不能只用 341 个差异蛋白
##   2)（可选）小鼠 TLE 蛋白质组全量表（文献数据），用于跨物种验证
##   3)（可选）ExoCarta/Vesiclepedia 外泌体蛋白列表文件（csv/txt，一列基因符号）
##      若未提供：脚本先尝试联网下载 ExoCarta，失败则用内置 Top100+ 列表
##
## 输出：results/23_exosome_proteome/ 下的 CSV + PDF
## =============================================================================
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root

dir.create(file.path(BASE_DIR, "ref"), showWarnings = FALSE)
for (f in c("ExoCarta_protein_genes_release6.csv", "ExoCarta_protein_mRNA_details_6.txt")) {
  file.copy(file.path(BASE_DIR, "BASE_DIR", "ref", f), file.path(BASE_DIR, "ref", f), overwrite = TRUE)
}
file.exists(file.path(BASE_DIR, "ref", "ExoCarta_protein_genes_release6.csv"))  # 必须 TRUE
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2)
  library(tibble); library(tidyr); library(stringr)
})

## --------------------------- 0. 配置区（按实际情况修改） ---------------------
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root

## 大鼠蛋白质组全量定量表（csv/xlsx 均可；xlsx 需 readxl）
RAT_PROTEOME_FILE <- file.path(BASE_DIR, "Differentially_expressed_statistics.xlsx")  # 大鼠KA蛋白组
RAT_SHEET         <- "KAvsnc"   # 该 xlsx 的 sheet 名
RAT_SPECIES       <- "rat"      # rat / mouse / human

## 小鼠 TLE 蛋白质组（文献补充材料；没有就留 NULL，脚本自动跳过跨物种部分）
MOUSE_PROTEOME_FILE <- file.path(BASE_DIR, "12974_2023_2840_MOESM3_ESM.xlsx")
MOUSE_SHEET         <- 1
MOUSE_HDR           <- 15      # 该文件前 14 行为图例说明，真正表头在第 15 行

## 外泌体蛋白列表（已下载好放在 ref/；留 NULL 则自动下载/用内置列表）
EXOCARTA_FILE <- file.path(BASE_DIR, "ref", "ExoCarta_protein_genes_release6.csv")

## 差异蛋白筛选阈值（与原分析保持一致）
PADJ_CUT  <- 0.05
LOG2FC_CUT <- log2(1.3)         # 若原分析用别的阈值，在这里改

OUT_DIR <- file.path(BASE_DIR, "results", "23_exosome_proteome")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## --------------------------- 1. 小工具函数 ----------------------------------
msg <- function(...) cat(sprintf(...), "\n")

## 防御性空表（避免 bind_rows 时 direction 变 logical(0) 的老坑）
empty_de <- function() tibble(gene = character(), log2FC = numeric(),
                              padj = numeric(), significant = logical())

## 读取表格（csv / tsv / xlsx 自适应；xlsx 可指定 sheet；hdr 可手动指定表头行号）
read_table_auto <- function(path, sheet = 1, hdr = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("请先 install.packages('readxl')")
    if (is.null(hdr)) {
      ## 智能表头探测：前 40 行里找"像表头"的行（≥2 个单元格命中列名关键词）
      raw <- readxl::read_excel(path, sheet = sheet, col_names = FALSE, n_max = 40)
      kw <- "^(gene|genes|gene name|gene symbol|protein|protein accession|accession|description|ratio|log2fc|logfc|fold change|p value|p\\.?value|adj|fdr|q ?value|peptides|unique peptides|intensity)"
      score <- apply(raw, 1, function(r)
        sum(grepl(kw, tolower(trimws(as.character(r))), perl = TRUE), na.rm = TRUE))
      hdr <- which(score >= 2)[1]
      if (is.na(hdr)) hdr <- 1
    }
    return(as_tibble(readxl::read_excel(path, sheet = sheet, skip = hdr - 1)))
  }
  if (ext %in% c("tsv", "txt")) return(readr::read_tsv(path, show_col_types = FALSE))
  readr::read_csv(path, show_col_types = FALSE)
}

## 从任意列名的蛋白质组表里提取 gene / log2FC / padj
standardize_proteome <- function(df, label = "") {
  cn <- names(df); cn_low <- tolower(cn)
  pick <- function(cands) {
    hit <- cn[which(cn_low %in% cands)[1]]
    if (is.na(hit)) hit <- cn[str_detect(cn_low, cands[1])][1]
    hit
  }
  g_col <- pick(c("gene", "gene name", "gene_name", "gene symbol", "gene_symbol",
                  "genesymbol", "symbol", "genename", "gene.names", "genes"))
  f_col <- pick(c("log2fc", "logfc", "log2foldchange", "log2 fold change", "log2(fc)"))
  p_col <- pick(c("p_adj", "padj", "adj.p.val", "adj p val", "fdr", "qvalue", "q value", "q_value"))
  if (is.na(p_col)) p_col <- pick(c("p value", "p_value", "pvalue", "p.val", "p"))
  if (is.na(p_col)) {  ## 兜底：任何含 "p value"/"p_value" 的列（如 "KA/nc P value"）
    p_col <- cn[str_detect(cn_low, "p.?value")][1]
  }
  ## log2FC 列缺失时，自动用 Ratio / Fold Change 列取 log2
  ratio_mode <- FALSE
  if (is.na(f_col)) {
    f_col <- cn[str_detect(cn_low, "ratio|fold.?change")][1]
    ratio_mode <- !is.na(f_col)
  }
  if (any(is.na(c(g_col, f_col, p_col))))
    stop(sprintf("[%s] 列名识别失败：gene=%s FC/Ratio=%s padj=%s；请把表头改成 gene/log2FC/padj",
                 label, g_col, f_col, p_col))
  if (ratio_mode) msg("[%s] 无 log2FC 列，自动对 '%s' 取 log2", label, f_col)
  out <- df %>%
    dplyr::transmute(gene   = as.character(.data[[g_col]]),
                     log2FC = if (ratio_mode) suppressWarnings(log2(as.numeric(.data[[f_col]])))
                     else as.numeric(.data[[f_col]]),
                     padj   = as.numeric(.data[[p_col]])) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(log2FC), is.finite(log2FC)) %>%
    dplyr::mutate(gene_upper = toupper(gene)) %>%
    dplyr::group_by(gene_upper) %>%
    dplyr::slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%   # 一基因多蛋白取最显著
    dplyr::ungroup()
  msg("[%s] 载入 %d 个可检测蛋白（去重后基因数）", label, nrow(out))
  out
}

mark_de <- function(df) {
  df %>% dplyr::mutate(padj_use = ifelse(is.na(padj), 1, padj),
                       significant = padj_use < PADJ_CUT & abs(log2FC) > LOG2FC_CUT)
}

## --------------------------- 2. 外泌体蛋白列表 ------------------------------
## 内置 ExoCarta Top100+ 经典外泌体蛋白（含 LGALS3 / CD9 / CD63 / CD81 / TSG101）
## 作为下载失败时的兜底；正式投稿建议换成 ExoCarta 全量列表
builtin_exosome <- c(
  "CD9","CD63","CD81","TSG101","PDCD6IP","SDCBP","HGS","STAM1","CHMP4A","CHMP4B",
  "VPS4A","VPS4B","VTA1","CHMP2A","RAB27A","RAB27B","RAB11A","RAB11B","RAB7A",
  "RAB5A","RAB1A","RAB2A","RAB8A","RAB10","RAB14","RAB35","ARF6","SMPD3","SMPD2",
  "HSPA8","HSPA1A","HSPA5","HSP90AA1","HSP90AB1","HSPB1","HSPD1","ENO1","ENO2",
  "GAPDH","PGK1","PKM","LDHA","LDHB","ALDOA","TPI1","GPI","MDH1","MDH2","CS",
  "ACTB","ACTG1","ACTA2","TUBA1A","TUBA1B","TUBB","TUBB4B","CFL1","PFN1","MSN",
  "EZR","RDX","ANXA1","ANXA2","ANXA5","ANXA6","ANXA7","ANXA11","S100A6","S100A8",
  "S100A9","S100A10","S100A11","S100A4","FLOT1","FLOT2","CAV1","STOM","BASP1",
  "GNAI2","GNB1","GNB2","GNA12","RHOA","RAC1","CDC42","ARF1","CLTC","AP2A1",
  "AP2B1","AP2M1","DNM2","ITGB1","ITGAV","ITGA2","ITGA5","ITGA6","ICAM1","VCAM1",
  "LGALS3","LGALS1","LGALS8","LGALS9","MVP","PRDX1","PRDX2","PRDX6","TXN","GSTP1",
  "PARK7","PGAM1","YWHAE","YWHAG","YWHAH","YWHAQ","YWHAZ","YWHAB","14-3-3",
  "C1QA","C1QB","C1QC","C3","C4A","SERPING1","APOE","APOA1","CLU","VTN","FGA",
  "FGB","FGG","ALB","TF","TTR","A2M","IGHG1","IGKC","HLA-A","HLA-B","B2M","CD47",
  "CD55","CD59","CD46","PTGFRN","BSG","SLC3A2","SLC1A5","SLC2A1","ATP1A1","ATP1B1",
  "VDAC1","VDAC2","UQCRC2","COX5B","ATP5A1","ATP5B","EEF1A1","EEF2","RPL3","RPS3",
  "HNRNPA2B1","HNRNPK","YBX1","SYNCRIP","AHNAK","KRT1","KRT10","KRT14","DSP"
)
builtin_exosome <- setdiff(unique(builtin_exosome), "14-3-3")  # 去掉占位项

get_exosome_list <- function() {
  if (!is.null(EXOCARTA_FILE) && file.exists(EXOCARTA_FILE)) {
    exo <- read_table_auto(EXOCARTA_FILE)
    cn_low <- tolower(names(exo))
    ## 自适应 ExoCarta/Vesiclepedia 官方下载格式：
    ##   含 "Content Type" 列时只保留 protein 行；基因列自动定位
    ct_idx <- which(str_detect(cn_low, "content.?type"))[1]
    if (!is.na(ct_idx)) exo <- exo %>% dplyr::filter(str_detect(tolower(as.character(.data[[ct_idx]])), "protein"))
    g_idx <- which(str_detect(cn_low, "gene"))[1]
    if (is.na(g_idx)) g_idx <- 1
    v <- toupper(unique(as.character(exo[[g_idx]])))
    v <- v[!is.na(v) & v != ""]
    msg("外泌体列表：来自文件 %s（%d 个基因）", basename(EXOCARTA_FILE), length(v))
    return(v)
  }
  dl <- tryCatch({
    if (!requireNamespace("httr", quietly = TRUE)) stop("no httr")
    ## ExoCarta 下载页可能变动；失败自动走内置列表
    tmp <- tempfile(fileext = ".csv")
    httr::GET("https://www.exocarta.org/exosome_markers/Archive/ExoCarta_protein_mRNA_details_6.txt",
              httr::timeout(60), httr::write_disk(tmp, overwrite = TRUE),
              httr::user_agent("Mozilla/5.0"))
    exo <- readr::read_csv(tmp, show_col_types = FALSE)
    cand <- names(exo)[tolower(names(exo)) %in% c("gene_symbol","genesymbol","symbol","entrez_gene","name")][1]
    if (is.na(cand)) cand <- names(exo)[2]
    v <- toupper(unique(as.character(exo[[cand]])))
    v <- v[!is.na(v) & v != ""]
    if (length(v) < 200) stop("downloaded list too small")
    v
  }, error = function(e) NULL)
  if (!is.null(dl)) {
    msg("外泌体列表：ExoCarta 在线下载成功（%d 个基因）", length(dl))
    return(dl)
  }
  msg("外泌体列表：下载失败，使用内置 Top100+ 列表（%d 个基因）；", length(builtin_exosome))
  msg("  正式投稿前建议到 exocarta.org / microvesicles.org 手动下载全量列表替换")
  builtin_exosome
}

## --------------------------- 3. 核心富集分析 -------------------------------
run_exosome_enrichment <- function(prot, exo, label, species = "rat") {
  prot <- mark_de(prot)

  ## ---- 降级模式：输入表只有差异蛋白（<1000 行）时，改用基因组近似背景 ----
  ## 注意：基因组背景会放大高丰度蛋白（外泌体 marker 普遍高丰度）的富集信号，
  ##       P 值只能作参考，论文中须声明此局限；找到全量定量表后应重跑标准模式。
  GENOME_BG <- c(rat = 22000, mouse = 21000, human = 20000)
  degraded <- nrow(prot) < 1000
  if (degraded) {
    msg("[%s] !! 输入仅 %d 行（差异蛋白表），启用降级模式：背景 = %s 基因组（%d 蛋白编码基因）",
        label, nrow(prot), species, GENOME_BG[[species]])
    bg  <- prot$gene_upper           # 仅占位
    de  <- unique(prot$gene_upper)   # 全表即差异蛋白，不再按阈值过滤
    exo_in_bg <- exo                 # 无法获知可检测性，用整个外泌体列表
    N <- GENOME_BG[[species]]
  } else {
    bg  <- unique(prot$gene_upper)                                   # 背景 = 全部可检测
    de  <- unique(prot$gene_upper[prot$significant])
    exo_in_bg <- intersect(exo, bg)                                  # 只统计可检测到的
    N <- length(bg)
  }
  hit   <- intersect(de, exo_in_bg)

  K <- length(exo_in_bg)
  n <- length(de); k <- length(hit)
  p_hyper <- if (k > 0) phyper(k - 1, K, N - K, n, lower.tail = FALSE) else 1
  msg("[%s] 背景 %d | 差异蛋白 %d | 外泌体蛋白(参照) %d | 命中 %d | 超几何 P = %.3g%s",
      label, N, n, K, k, p_hyper, ifelse(degraded, "（降级模式，仅供参考）", ""))

  res <- list(
    summary = tibble(dataset = label, mode = ifelse(degraded, "degraded_genome_bg", "standard_detected_bg"),
                     background = N, DE = n, exosome_detected = K,
                     overlap = k, fold_enrichment = round((k / max(n,1)) / (K / N), 2),
                     hypergeometric_p = signif(p_hyper, 3)),
    hit_genes = prot %>% dplyr::filter(gene_upper %in% hit) %>%
      dplyr::arrange(padj_use) %>% dplyr::select(gene, log2FC, padj),
    prot = prot, exo_in_bg = exo_in_bg, hit = hit
  )

  ## 关键分子检查
  key <- c("LGALS3","CD9","CD63","CD81","TSG101","PDCD6IP","SDCBP","TREM2")
  key_tab <- tibble(gene = key,
                    in_exosome_list = key %in% exo,
                    detected_in_proteome = key %in% bg,
                    DE = key %in% de,
                    log2FC = sapply(key, function(g) {
                      v <- prot$log2FC[prot$gene_upper == g]; if (length(v)) round(v[1],3) else NA
                    }))
  msg("[%s] 关键分子：", label); print(as.data.frame(key_tab))
  res$key <- key_tab

  ## 糖基化/凝集素相关外泌体蛋白
  glyco <- intersect(exo_in_bg, c("LGALS1","LGALS3","LGALS8","LGALS9","C1QA","C1QB","C1QC",
                                  "MAN1A1","MAN2A1","MGAT1","MGAT5","B4GALT1","ST6GAL1",
                                  "ST3GAL4","FUT8","GANAB","GLG1","CALR","CANX"))
  res$glyco <- prot %>% dplyr::filter(gene_upper %in% glyco) %>%
    dplyr::select(gene, log2FC, padj) %>% dplyr::arrange(padj)

  res
}

## --------------------------- 4. 绘图 ---------------------------------------
plot_venn <- function(res, label) {
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    msg("未安装 VennDiagram，跳过韦恩图（install.packages('VennDiagram')）"); return(invisible())
  }
  de <- unique(res$prot$gene_upper[res$prot$significant])
  vd <- VennDiagram::venn.diagram(
    x = list(`Differential proteins` = de, `Exosome proteins (detected)` = res$exo_in_bg),
    filename = NULL, fill = c("#C9828B", "#2E86C8"), alpha = 0.55,
    cex = 1.4, cat.cex = 1.2, cat.fontface = "bold", margin = 0.08)
  pdf(file.path(OUT_DIR, sprintf("23_venn_%s.pdf", label)), width = 5.5, height = 5)
  grid::grid.draw(vd); grid::grid.text(label, y = 0.97, gp = grid::gpar(fontface = "bold"))
  dev.off()
  msg("[%s] 韦恩图已保存", label)
}

plot_heatmap <- function(res, label, topn = 40) {
  df <- res$hit_genes %>% dplyr::slice_head(n = topn)
  if (nrow(df) < 3) { msg("[%s] 命中过少，跳过热图", label); return(invisible()) }
  p <- ggplot(df %>% dplyr::mutate(gene = factor(gene, levels = rev(gene))),
              aes(x = "KA vs Con", y = gene, fill = log2FC)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", log2FC)), size = 3) +
    scale_fill_gradient2(low = "#2E86C8", mid = "white", high = "#C9828B", midpoint = 0) +
    labs(title = sprintf("Exosome proteins among DE proteins — %s (top %d)", label, nrow(df)),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11),
          axis.text.y = element_text(size = 8))
  ggsave(file.path(OUT_DIR, sprintf("23_heatmap_%s.pdf", label)), p,
         width = 4.5, height = max(3.5, 0.22 * nrow(df) + 1.5))
  msg("[%s] 热图已保存（%d 个蛋白）", label, nrow(df))
}

plot_go <- function(res, label, species) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    msg("未安装 clusterProfiler，跳过 GO/KEGG"); return(invisible())
  }
  orgdb <- switch(species,
                  rat   = "org.Rn.eg.db", mouse = "org.Mm.eg.db", human = "org.Hs.eg.db")
  if (!requireNamespace(orgdb, quietly = TRUE)) {
    msg("未安装 %s，跳过 GO/KEGG", orgdb); return(invisible())
  }
  OrgDb <- getExportedValue(orgdb, orgdb)
  ## 标准模式用"可检测蛋白"做 universe 避免高丰度偏倚；降级模式（仅差异表）不设 universe
  uni <- res$prot %>% dplyr::filter(!is.na(log2FC))
  use_universe <- nrow(uni) >= 1000
  bg_id <- if (use_universe) clusterProfiler::bitr(uni$gene_upper, fromType = "SYMBOL",
                                                   toType = "ENTREZID", OrgDb = OrgDb) else NULL
  hit_id <- clusterProfiler::bitr(res$hit, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  if (nrow(hit_id) < 5) { msg("[%s] 命中基因映射后 <5，跳过富集", label); return(invisible()) }
  ego <- if (use_universe) {
    clusterProfiler::enrichGO(hit_id$ENTREZID, OrgDb = OrgDb, ont = "BP",
                              universe = bg_id$ENTREZID,
                              pAdjustMethod = "BH", qvalueCutoff = 0.2)
  } else {
    clusterProfiler::enrichGO(hit_id$ENTREZID, OrgDb = OrgDb, ont = "BP",
                              pAdjustMethod = "BH", qvalueCutoff = 0.2)
  }
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write_csv(as.data.frame(ego), file.path(OUT_DIR, sprintf("23_GO_%s.csv", label)))
    p <- clusterProfiler::dotplot(ego, showCategory = 15) +
      ggtitle(sprintf("GO BP of exosome DE proteins — %s", label))
    ggsave(file.path(OUT_DIR, sprintf("23_GO_%s.pdf", label)), p, width = 8, height = 6)
    msg("[%s] GO 富集已保存（%d 条通路）", label, nrow(as.data.frame(ego)))
  } else msg("[%s] GO 无显著富集", label)
}

## --------------------------- 5. 主流程 -------------------------------------
exo_list <- get_exosome_list()
write_csv(tibble(exosome_gene = exo_list), file.path(OUT_DIR, "23_exosome_list_used.csv"))

rat <- standardize_proteome(read_table_auto(RAT_PROTEOME_FILE, sheet = RAT_SHEET), label = "rat_KA")
res_rat <- run_exosome_enrichment(rat, exo_list, "rat_KA", species = RAT_SPECIES)
write_csv(res_rat$summary,   file.path(OUT_DIR, "23_enrichment_summary_rat.csv"))
write_csv(res_rat$hit_genes, file.path(OUT_DIR, "23_exosome_DE_hits_rat.csv"))
write_csv(res_rat$key,       file.path(OUT_DIR, "23_key_molecules_rat.csv"))
write_csv(res_rat$glyco,     file.path(OUT_DIR, "23_glyco_exosome_rat.csv"))
plot_venn(res_rat, "rat_KA"); plot_heatmap(res_rat, "rat_KA"); plot_go(res_rat, "rat_KA", RAT_SPECIES)

## --------------------------- 6. 跨物种验证（可选） --------------------------
if (!is.null(MOUSE_PROTEOME_FILE) && file.exists(MOUSE_PROTEOME_FILE)) {
  mou <- standardize_proteome(read_table_auto(MOUSE_PROTEOME_FILE, sheet = MOUSE_SHEET, hdr = MOUSE_HDR), label = "mouse_TLE")
  res_mou <- run_exosome_enrichment(mou, exo_list, "mouse_TLE", species = "mouse")
  write_csv(res_mou$summary,   file.path(OUT_DIR, "23_enrichment_summary_mouse.csv"))
  write_csv(res_mou$hit_genes, file.path(OUT_DIR, "23_exosome_DE_hits_mouse.csv"))
  write_csv(res_mou$key,       file.path(OUT_DIR, "23_key_molecules_mouse.csv"))
  plot_venn(res_mou, "mouse_TLE"); plot_heatmap(res_mou, "mouse_TLE")
  plot_go(res_mou, "mouse_TLE", "mouse")

  shared <- intersect(res_rat$hit, res_mou$hit)
  msg("跨物种共享外泌体差异蛋白 %d 个：%s", length(shared), paste(shared, collapse = ", "))
  write_csv(tibble(shared_exosome_DE = shared),
            file.path(OUT_DIR, "23_cross_species_shared.csv"))
  write_csv(bind_rows(res_rat$summary, res_mou$summary),
            file.path(OUT_DIR, "23_enrichment_summary_all.csv"))

  if (requireNamespace("VennDiagram", quietly = TRUE)) {
    vd <- VennDiagram::venn.diagram(
      x = list(`Rat KA exosome DE` = res_rat$hit, `Mouse TLE exosome DE` = res_mou$hit),
      filename = NULL, fill = c("#684765", "#D9A441"), alpha = 0.55,
      cex = 1.4, cat.cex = 1.2, cat.fontface = "bold", margin = 0.08)
    pdf(file.path(OUT_DIR, "23_venn_cross_species.pdf"), width = 5.5, height = 5)
    grid::grid.draw(vd); grid::grid.text("Cross-species conservation",
                                         y = 0.97, gp = grid::gpar(fontface = "bold"))
    dev.off()
  }
} else {
  msg("未提供小鼠蛋白质组表，跳过跨物种验证（MOUSE_PROTEOME_FILE 留空）")
}

msg("===== 全部完成，输出目录：%s =====", OUT_DIR)
