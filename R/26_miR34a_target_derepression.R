# ============================================================================
# Script      : 26_miR34a_target_derepression.R
# Purpose     : Collective de-repression test of miR-34a-5p targets in TLE mRNA cohorts (pathway-level support for miRNA down -> target up).
# Manuscript  : Supplementary analyses (miRNA axis)
# Input       : GSE168375 mRNA logFC; multiMiR target list
# Output      : results/26_* test tables/figures
# Dependencies: dplyr, ggplot2, readr, stringr, tibble
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## =============================================================================
## 26_miR34a_target_derepression.R
## 方案 3.2 加分项：miR-34a-5p 靶基因在 TLE mRNA 队列中是否整体去抑制（上调）？
##
## 逻辑：
##   miR-34a-5p 在 TLE 海马下调（脚本 25，GSE99455 FDR=0.018，GSE245228 同向）
##   → 若其靶基因整体上调，则为「miRNA 下调 → 靶 mRNA 去抑制」提供通路级证据
##   → LGALS3 只是其中一个靶基因，该分析把单基因观察升级为机制链条
##
## 数据：
##   mRNA 队列：GSE168375（人皮层 bulk，31 TLE vs 12 Con）——与脚本 24 同一文件
##   靶基因清单：multiMiR（miRTarBase/TarBase 实验验证 + ≥2 预测库）
##
## 方法：
##   1) limma/voom 或 Wilcoxon 算全部基因 TLE vs Con 的 logFC
##   2) 靶基因 vs 非靶基因的 logFC 分布比较（Wilcoxon 秩和，单尾：靶 > 非靶）
##   3) 置换检验（1000 次随机同大小基因集）得经验 P 值
##   4) 输出：靶基因明细表、LGALS3 在其中的位置、分布图
##
## 输出：results/26_miR34a_derepression/
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(stringr); library(ggplot2)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR  <- file.path(BASE_DIR, "results", "26_miR34a_derepression")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(sprintf(...), "\n")

## 【配置】与脚本 24 相同的 GSE168375 表达矩阵文件
GSE168375_FILE <- file.path(DATA_DIR, "GSE168375_expr.csv")  # 与脚本 24 同一文件（已由 24 自动生成）

## --------------------------- 1. miR-34a-5p 靶基因 ----------------------------
get_mir34a_targets <- function() {
  if (requireNamespace("multiMiR", quietly = TRUE)) {
    q <- tryCatch(suppressMessages(
      multiMiR::get_multimir(org = "hsa", mirna = "hsa-miR-34a-5p",
                             table = "all", summary = TRUE)),
      error = function(e) { msg("multiMiR 查询报错：%s", conditionMessage(e)); NULL })
    if (!is.null(q)) {
      d <- q@data
      msg("multiMiR 原始返回 %d 行；miRNA 分布前 5：%s", nrow(d),
          paste(head(names(sort(table(d$mature_mirna_id), decreasing = TRUE)), 5), collapse = ", "))
      ## 手动按 miRNA 名过滤（mirna 参数本身过滤不可靠，双保险）
      d <- d %>% dplyr::filter(mature_mirna_id %in% c("hsa-miR-34a-5p","hsa-miR-34a"))
      msg("过滤后 miR-34a-5p 互作 %d 条，涉及库：%s", nrow(d),
          paste(unique(d$database), collapse = ", "))
      ## 实验验证（miRTarBase/TarBase，优先强证据）
      val <- d %>% dplyr::filter(database %in% c("mirtarbase","tarbase")) %>%
        dplyr::pull(target_symbol) %>% unique()
      ## 预测库 ≥2 个支持
      pred <- d %>% dplyr::filter(database %in% c("mirdb","targetscan","diana_microt")) %>%
        dplyr::count(target_symbol, name = "n_pred") %>%
        dplyr::filter(n_pred >= 2) %>% dplyr::pull(target_symbol)
      tg <- union(val, pred)
      tg <- tg[!is.na(tg) & tg != ""]
      if (length(tg) >= 20) {
        msg("miR-34a-5p 靶基因 %d 个（实验验证 %d 个 + 预测≥2库 %d 个）",
            length(tg), length(val), length(pred))
        return(list(all = tg, validated = val))
      }
    }
  }
  msg("!! multiMiR 失败，走 26 基因兜底清单——结果仅供探索，不作结论依据")
  ## 经典验证靶（文献高频报道）
  list(all = c("LGALS3","SIRT1","BCL2","MYC","CDK4","CDK6","CCND1","CCNE2","E2F1","E2F3",
               "MET","AXL","PDGFRB","NOTCH1","NOTCH2","DLL1","JAG1","HMGA2","FOXP1",
               "TNFSF11","WNT1","WNT3","SATB2","GRM7","CACNA1C","SYT1"),
       validated = character(0))
}

## 【默认】TargetScan 官方离线清单（绕开 multiMiR——其查询结果两次验证不可靠）
## 文件由 TargetScan 8.0 Predicted_Targets_Info 提取（miR-34-5p/449-5p 家族，人源）
TS_FILE <- file.path(BASE_DIR, "ref", "TargetScan_miR-34a-5p_targets.csv")
if (file.exists(TS_FILE)) {
  ts <- read_csv(TS_FILE, show_col_types = FALSE)
  msg("TargetScan miR-34-5p/449-5p 家族靶基因 %d 个（高置信 %d 个）",
      nrow(ts), sum(ts$high_confidence))
  tg <- list(all = ts$`Gene Symbol`,
             validated = ts$`Gene Symbol`[ts$high_confidence])
} else {
  msg("未找到 %s，回退 multiMiR/兜底清单", TS_FILE)
  tg <- get_mir34a_targets()
}
write_csv(tibble(target_gene = tg$all, experimentally_validated = tg$all %in% tg$validated),
          file.path(OUT_DIR, "26_miR34a5p_targets.csv"))

## --------------------------- 2. 读取 mRNA 表达矩阵 ---------------------------
if (!file.exists(GSE168375_FILE)) {
  stop("找不到 GSE168375 文件：", GSE168375_FILE,
       "\n请修改脚本顶部 GSE168375_FILE 指向你的表达矩阵（与脚本 24 同一文件）")
}

read_expr <- function(path) {
  df <- if (str_detect(path, "\\.(csv|txt|tsv)(\\.gz)?$")) {
    readr::read_delim(path, delim = NULL, show_col_types = FALSE)
  } else readxl::read_excel(path)
  msg("表达矩阵：%d 行 x %d 列；前 5 列：%s", nrow(df), ncol(df),
      paste(head(names(df), 5), collapse = ", "))
  idcol <- names(df)[1]
  mat <- df %>%
    dplyr::mutate(gene = as.character(.data[[idcol]])) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  m <- as.matrix(mat[, -1, drop = FALSE]); rownames(m) <- mat$gene
  m
}

expr <- read_expr(GSE168375_FILE)

## --------------------------- 3. 分组与差异 -----------------------------------
## GSE168375：31 TLE vs 12 Con。列名为 GSM 编号，须用 GEO 官方映射分组
## （gsm_map 在脚本 24 中定义；若当前会话没有则重建）
if (!exists("gsm_map")) {
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
}
dx <- gsm_map$dx[match(colnames(expr), gsm_map$gsm)]
grp <- ifelse(dx == "epilepsy", "TLE", "Con"); names(grp) <- colnames(expr)
msg("分组：TLE=%d Con=%d（应为 31 vs 12）", sum(grp=="TLE"), sum(grp=="Con"))
stopifnot(sum(grp=="TLE") == 31, sum(grp=="Con") == 12)

## 低表达过滤（在 log 转换前，用 counts 尺度）：至少 10 个样本 CPM >= 1
cpm <- t(t(expr) / colSums(expr) * 1e6)
keep_g <- rowSums(cpm >= 1) >= 10
expr <- expr[keep_g, ]
msg("低表达过滤：%d / %d 基因保留", sum(keep_g), length(keep_g))

## log2 转换检测（最大值 > 100 认为未取 log）
if (max(expr, na.rm = TRUE) > 100) {
  expr <- log2(expr + 1)
  msg("已执行 log2(x+1) 转换")
}

## 全基因 logFC（TLE 均值 - Con 均值）+ Wilcoxon P
gene_de <- lapply(rownames(expr), function(g) {
  x <- expr[g, ]
  p <- tryCatch(wilcox.test(x[grp=="TLE"], x[grp=="Con"])$p.value, error = function(e) NA)
  tibble(gene = g, logFC = mean(x[grp=="TLE"]) - mean(x[grp=="Con"]), pval = p)
}) %>% bind_rows() %>%
  dplyr::mutate(fdr = p.adjust(pval, "BH"),
                is_target = gene %in% tg$all,
                is_validated_target = gene %in% tg$validated)

write_csv(gene_de, file.path(OUT_DIR, "26_gene_de_with_target_flag.csv"))

## --------------------------- 4. 去抑制检验 -----------------------------------
de_t  <- gene_de %>% dplyr::filter(is_target, is.finite(logFC))
de_nt <- gene_de %>% dplyr::filter(!is_target, is.finite(logFC))
msg("靶基因 %d 个在矩阵中检出，非靶基因 %d 个", nrow(de_t), nrow(de_nt))

## 4a. Wilcoxon 秩和（单尾：靶基因 logFC 更高 = 上调）
wt <- wilcox.test(de_t$logFC, de_nt$logFC, alternative = "greater")
msg("Wilcoxon（靶 > 非靶，单尾）：P = %.4g；靶基因中位 logFC = %.3f，非靶 = %.3f",
    wt$p.value, median(de_t$logFC), median(de_nt$logFC))

## 4b. 置换检验（经验 P）
set.seed(42)
obs  <- median(de_t$logFC)
perm <- replicate(1000, median(sample(de_nt$logFC, nrow(de_t))))
emp_p <- (sum(perm >= obs) + 1) / 1001
msg("置换检验（1000 次）：经验 P = %.4f（观测中位 logFC = %.3f，随机中位 = %.3f）",
    emp_p, obs, median(perm))

## 4c. LGALS3 位置
lg <- gene_de %>% dplyr::filter(gene == "LGALS3")
if (nrow(lg) > 0) {
  pct <- mean(de_nt$logFC <= lg$logFC[1]) * 100
  msg("LGALS3：logFC = %.3f，P = %.4g，FDR = %.4g（高于 %.1f%% 的非靶基因）",
      lg$logFC, lg$pval, lg$fdr, pct)
}

## 4d. 实验验证靶基因单独再验一次（强证据子集）
if (length(tg$validated) >= 10) {
  de_v <- gene_de %>% dplyr::filter(is_validated_target, is.finite(logFC))
  wt_v <- wilcox.test(de_v$logFC, de_nt$logFC, alternative = "greater")
  msg("[验证靶子集 %d 个] Wilcoxon 单尾 P = %.4g，中位 logFC = %.3f",
      nrow(de_v), wt_v$p.value, median(de_v$logFC))
}

## --------------------------- 5. 判定与图 -------------------------------------
verdict <- if (wt$p.value < 0.05 && emp_p < 0.05 && obs > 0) {
  ">>> 判定：支持去抑制——miR-34a-5p 靶基因在 TLE 中整体上调，通路级证据成立。\n>>> 写法：miR-34a-5p 下调 → 靶基因（含 LGALS3）整体去抑制 → 部分解释 LGALS3 蛋白升高。"
} else if (obs > 0 && wt$p.value < 0.1) {
  ">>> 判定：弱支持——靶基因有上调趋势但未达显著，写为提示性证据，勿下强结论。"
} else {
  ">>> 判定：不支持——靶基因无整体上调，miR-34a-5p 下调不足以解释 LGALS3，转向翻译/降解/外泌体层面。"
}
cat(verdict, "\n")

plot_df <- bind_rows(
  de_nt %>% dplyr::sample_n(min(5000, n())) %>% dplyr::mutate(set = "Non-targets"),
  de_t %>% dplyr::mutate(set = "miR-34a-5p targets")
)
p <- ggplot(plot_df, aes(x = set, y = logFC, fill = set)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  {if (nrow(lg) > 0) annotate("point", x = 2, y = lg$logFC[1], color = "#D9A441", size = 3)} +
  {if (nrow(lg) > 0) annotate("text", x = 2, y = lg$logFC[1], label = "LGALS3",
                              hjust = -0.2, color = "#D9A441", fontface = "bold")} +
  scale_fill_manual(values = c("Non-targets" = "grey70", "miR-34a-5p targets" = "#0E3A6B")) +
  labs(title = "De-repression of miR-34a-5p targets in TLE",
       subtitle = sprintf("Wilcoxon one-sided P = %.3g; permutation P = %.3f (n=%d targets)",
                          wt$p.value, emp_p, nrow(de_t)),
       y = "logFC (TLE vs Con, mRNA)", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")
ggsave(file.path(OUT_DIR, "26_miR34a_derepression.pdf"), p, width = 5.5, height = 5)
msg("图已保存：%s", file.path(OUT_DIR, "26_miR34a_derepression.pdf"))

msg("===== 完成，输出目录：%s =====", OUT_DIR)
