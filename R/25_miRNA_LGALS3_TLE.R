# ============================================================================
# Script      : 25_miRNA_LGALS3_TLE.R
# Purpose     : miRNA regulation: are LGALS3-targeting miRNAs down-regulated in TLE hippocampus (GSE99455 main, GSE245228 validation)?
# Manuscript  : Supplementary analyses (miRNA axis)
# Input       : GSE99455 miRNA-seq counts; GSE245228; multiMiR (online) or manual target list
# Output      : results/25_* candidate miRNA-LGALS3 axis tables
# Dependencies: dplyr, ggplot2, readr, stringr, tibble
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# Setup       : package-install calls commented out ([setup] tags); run 00_setup_genesets.R first
# ============================================================================

## =============================================================================
## 25_miRNA_LGALS3_TLE.R
## 方案 3.2：miRNA 调控——靶向 LGALS3 的 miRNA 在 TLE 海马中是否下调？
##
## 数据：
##   主队列  GSE99455（人海马 miRNA-seq，16 mTLE+HS vs 8 对照）
##   验证队列 GSE245228 海马胞质组分（5v5，注意其分组为 CTL vs W0）
##
## 流程：
##   1) 读 counts → 差异 miRNA（edgeR；未安装则 Wilcoxon on log-CPM 兜底）
##   2) LGALS3 靶向 miRNA 清单：multiMiR 在线查询（miRDB/TargetScan/miRTarBase），
##      仅保留 ≥2 库支持或有实验验证的；无网络/无包时用手工清单（需你确认）
##   3) 交集：下调 miRNA ∩ LGALS3 靶向 miRNA → 候选调控轴
##   4) （加分）miRNA 靶基因去抑制：候选 miRNA 的全部靶基因在 mRNA 队列中是否上调
##
## 输出：results/25_miRNA/
## =============================================================================
# [setup] one-time package installation commented out for submission;
# [setup] run 00_setup_genesets.R (or install manually) before first use:
# BiocManager::install(c("edgeR", "multiMiR"))
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(stringr); library(ggplot2)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "results", "25_miRNA")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(sprintf(...), "\n")

## --------------------------- 1. 定位输入文件 --------------------------------
find1 <- function(pattern) {
  f <- list.files(DATA_DIR, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) NULL else f[1]
}



find1 <- function(pattern, exclude = NULL) {
  f <- list.files(DATA_DIR, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (!is.null(exclude)) f <- f[!grepl(exclude, basename(f), ignore.case = TRUE)]
  if (length(f) == 0) NULL else f[1]
}

F_99455  <- find1("GSE99455.*plain_counts_per_file", exclude = "Normali")
F_99455N <- find1("GSE99455.*Normali")
F_245cyt <- find1("GSE245228_Hip_cyt_counttable", exclude = "isomiRs")

cat("99455 原始counts :", basename(F_99455), "\n")
cat("99455 标准化表   :", basename(F_99455N), "\n")
cat("245228 海马胞质  :", basename(F_245cyt), "\n")
## --------------------------- 2. 读取 counts ---------------------------------
read_counts <- function(path, label) {
  df <- readr::read_tsv(path, show_col_types = FALSE)
  msg("[%s] 原始表 %d 行 x %d 列；前 6 列名：%s", label, nrow(df), ncol(df),
      paste(head(names(df), 6), collapse = ", "))
  idcol <- names(df)[1]
  mat <- df %>%
    dplyr::mutate(mir = as.character(.data[[idcol]])) %>%
    dplyr::filter(!is.na(mir), mir != "") %>%
    dplyr::group_by(mir) %>%
    dplyr::summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  m <- as.matrix(mat[, -1, drop = FALSE]); rownames(m) <- mat$mir
  msg("[%s] %d 个 miRNA x %d 个样本", label, nrow(m), ncol(m))
  m
}

## --------------------------- 3. 差异分析 ------------------------------------
de_mirna <- function(mat, grp, label) {
  grp <- grp[colnames(mat)]
  if (requireNamespace("edgeR", quietly = TRUE)) {
    y <- edgeR::DGEList(counts = mat, group = grp)
    keep <- edgeR::filterByExpr(y)
    y <- edgeR::calcNormFactors(y[keep, ])
    des <- model.matrix(~grp)
    y <- edgeR::estimateDisp(y, des)
    fit <- edgeR::glmQLFit(y, des); qlf <- edgeR::glmQLFTest(fit)
    tt <- edgeR::topTags(qlf, n = Inf)$table
    out <- tibble(mirna = rownames(tt), logFC = tt$logFC, pval = tt$PValue, fdr = tt$FDR,
                  method = "edgeR-QLF")
    msg("[%s] edgeR：%d 个 miRNA 入检，FDR<0.05 且下调 %d 个", label, nrow(out),
        sum(out$fdr < 0.05 & out$logFC < 0))
  } else {
    msg("[%s] 未安装 edgeR，用 Wilcoxon on log2-CPM 兜底（建议 BiocManager::install('edgeR')）", label)
    cpm <- t(t(mat) / colSums(mat) * 1e6); lcpm <- log2(cpm + 1)
    keep <- rowSums(cpm >= 1) >= 3; lcpm <- lcpm[keep, ]
    out <- lapply(rownames(lcpm), function(g) {
      x <- lcpm[g, ]
      p <- tryCatch(wilcox.test(x[grp == "TLE"], x[grp == "Con"])$p.value, error = function(e) NA)
      tibble(mirna = g, logFC = mean(x[grp == "TLE"]) - mean(x[grp == "Con"]), pval = p)
    }) %>% bind_rows() %>% dplyr::mutate(fdr = p.adjust(pval, "BH"), method = "wilcoxon-logCPM")
    msg("[%s] Wilcoxon：%d 个 miRNA 入检，FDR<0.05 且下调 %d 个", label, nrow(out),
        sum(out$fdr < 0.05 & out$logFC < 0, na.rm = TRUE))
  }
  out$dataset <- label
  out
}

## --------------------------- 4. LGALS3 靶向 miRNA ---------------------------
get_lgals3_mirs <- function() {
  if (requireNamespace("multiMiR", quietly = TRUE)) {
    q <- tryCatch({
      suppressMessages(multiMiR::get_multimir(org = "hsa", target = "LGALS3",
                                              table = "all", summary = TRUE))
    }, error = function(e) NULL)
    if (!is.null(q)) {
      d <- q@data
      d <- d %>% dplyr::filter(database %in% c("mirtarbase","tarbase") | .data[["support_type"]] != "")
      ## 保留：实验验证库命中 或 ≥2 个预测库命中
      pred <- d %>% dplyr::filter(database %in% c("mirdb","targetscan","diana_microt")) %>%
        dplyr::count(mature_mirna_acc, name = "n_pred")
      val  <- d %>% dplyr::filter(database %in% c("mirtarbase","tarbase"))
      keep_acc <- union(pred$mature_mirna_acc[pred$n_pred >= 2], unique(val$mature_mirna_acc))
      idmap <- d %>% dplyr::distinct(mature_mirna_acc, mature_mirna_id)
      v <- unique(idmap$mature_mirna_id[idmap$mature_mirna_acc %in% keep_acc])
      v <- v[!is.na(v)]
      if (length(v) >= 3) {
        msg("multiMiR：LGALS3 靶向 miRNA %d 个（≥2 库或实验验证）", length(v))
        return(v)
      }
    }
  }
  msg("multiMiR 不可用或结果过少，使用内置文献候选清单（请在 multiMiR 装好后复核）")
  msg("  安装：BiocManager::install('multiMiR')")
  ## 内置候选：文献报道与 LGALS3 3'UTR 相关的 miRNA（保守起见仅作探索）
  c("hsa-miR-21-5p","hsa-miR-22-3p","hsa-miR-128-3p","hsa-miR-146a-5p",
    "hsa-miR-155-5p","hsa-miR-199a-3p","hsa-miR-320a","hsa-miR-27a-3p")
}

## --------------------------- 5. 主流程 --------------------------------------
de_list <- list()

## -- GSE99455：GEO 元数据确认 S1-S16 = mTLE+HS，S17-S24 = 对照
if (!is.null(F_99455)) {
  m1 <- read_counts(F_99455, "GSE99455")
  snum <- as.numeric(str_extract(colnames(m1), "\\d+"))
  g1 <- ifelse(snum <= 16, "TLE", "Con"); names(g1) <- colnames(m1)
  msg("GSE99455 分组：TLE=%d Con=%d（应为 16 vs 8）", sum(g1=="TLE"), sum(g1=="Con"))
  if (sum(g1=="TLE") >= 5 && sum(g1=="Con") >= 4)
    de_list[["GSE99455"]] <- de_mirna(m1, g1, "GSE99455")
}

## -- GSE245228 海马胞质 miRNA counts：列顺序 CTLx5 然后 W0x5
F_245mir <- find1("GSE245228_Hip_cyt_CTLvsW0_countTable_edgeR", exclude = "isomiRs")
if (!is.null(F_245mir)) {
  m2 <- read_counts(F_245mir, "GSE245228_cyt")
  g2 <- rep(c("Con","TLE"), each = 5); names(g2) <- colnames(m2)
  msg("GSE245228 分组：TLE(W0)=%d Con(CTL)=%d", sum(g2=="TLE"), sum(g2=="Con"))
  de_list[["GSE245228_cyt"]] <- de_mirna(m2, g2, "GSE245228_cyt")
} else {
  msg("未找到 Hip_cyt_CTLvsW0_countTable_edgeR 文件，跳过 GSE245228（请先下载）")
}

if (length(de_list) == 0) stop("没有可用队列，请检查 data/ 下的文件")

de_all <- bind_rows(de_list)
write_csv(de_all, file.path(OUT_DIR, "25_DE_miRNA_all.csv"))

## -- 交集 --
targets <- get_lgals3_mirs()
write_csv(tibble(lgals3_targeting_mirna = targets),
          file.path(OUT_DIR, "25_LGALS3_targeting_miRNAs.csv"))

norm_mir <- function(x) {
  x <- tolower(x)
  x <- str_replace(x, "^hsa-", "")
  x
}
tg_norm <- norm_mir(targets)

hits <- de_all %>%
  dplyr::mutate(mir_norm = norm_mir(mirna),
                hit = mir_norm %in% tg_norm |
                      str_replace(mir_norm, "-[35]p$","") %in% str_replace(tg_norm, "-[35]p$","")) %>%
  dplyr::filter(hit) %>%
  dplyr::arrange(dataset, fdr)
write_csv(hits, file.path(OUT_DIR, "25_LGALS3_axis_hits.csv"))

msg("===== 靶向 LGALS3 的 miRNA 在 TLE 中的变化 =====")

if (nrow(hits) == 0) {
  msg("（无交集命中——若使用内置清单，建议先装 multiMiR 换完整预测清单重跑）")
} else {
  print(as.data.frame(hits %>% dplyr::select(dataset, mirna, logFC, pval, fdr)))

  down <- hits %>% dplyr::filter(logFC < 0, fdr < 0.1)
  if (nrow(down) > 0) {
    msg(">>> 核心候选（下调 + FDR<0.1）：%s", paste(down$mirna, collapse = ", "))
    msg(">>> 机制解释：miRNA 下调 → LGALS3 mRNA 去抑制 → 翻译增加/降解减少")
  } else {
    msg("交集中无显著下调者——miRNA 层面对 LGALS3 的解释力有限，转向翻译/降解/外泌体层面")
  }

  ## 先建图，再保存（顺序不能反）
  p <- ggplot(hits, aes(x = mirna, y = logFC, fill = dataset)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(yintercept = 0, color = "grey60") +
    coord_flip() +
    scale_fill_manual(values = c("GSE99455" = "#0E3A6B", "GSE245228_cyt" = "#D9A441")) +
    labs(title = "LGALS3-targeting miRNAs in TLE hippocampus",
         subtitle = "Negative logFC = down-regulated in TLE (supports de-repression of LGALS3)",
         y = "logFC (TLE vs Con)", x = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(OUT_DIR, "25_LGALS3_miRNA_hits.pdf"), p,
         width = 7.5, height = 49, limitsize = FALSE)
  msg("图已保存：%s", file.path(OUT_DIR, "25_LGALS3_miRNA_hits.pdf"))
}
msg("===== 完成，输出目录：%s =====", OUT_DIR)
