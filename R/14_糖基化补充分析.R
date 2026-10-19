# ============================================================================
# Script      : 14_糖基化补充分析.R
# Purpose     : Supplementary glycosylation analyses: lectin landscape, LGALS3 co-expression/bridge tests, glycosylation-enzyme panel, CellChat, blood.
# Manuscript  : Supplementary analyses (Galectin-3 axis related)
# Input       : existing results/ files only (no new downloads)
# Output      : results/14_* tables; glycosylation figures
# Dependencies: CellChat, GSVA, UCell, data.table, dplyr, ggplot2, pheatmap, readxl, tidyr
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ============================================================
## 14_糖基化补充分析.R（完整修正版）
## 凝集素全景 + LGALS3共表达/桥梁检验 + 糖基化酶谱 + CellChat + 血液
## 数据：全部使用已有文件，无需下载新数据集
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr)
  library(tidyr); library(pheatmap); library(readxl)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------- 0. 通用工具 ----------
find1 <- function(pattern){
  f <- list.files(OUT_DIR, pattern, full.names = TRUE)
  cat(sprintf("  [匹配 %s] -> %s\n", pattern, ifelse(length(f)>0, basename(f[1]), "【未找到】")))
  if(length(f)==0) return(NULL)
  f[1]
}

read_mat <- function(path){
  d <- fread(path, data.table = FALSE)
  m <- as.matrix(d[, -1, drop = FALSE])
  rownames(m) <- make.unique(as.character(d[[1]]))
  m
}

## 统一分组规则：CHKJ-H=疾病/CHKJ-C=对照；Drug-nave=疾病/Healthy=对照
split_grp <- function(cn){
  cn2 <- gsub("\\.", "-", cn)
  if(all(grepl("^CHKJ", cn2))){
    return(list(dis = cn[grepl("^CHKJ-H", cn2)], con = cn[grepl("^CHKJ-C", cn2)]))
  }
  if(any(grepl("Healthy", cn2))){
    return(list(dis = cn[grepl("Drug", cn2, ignore.case = TRUE)],
                con = cn[grepl("Healthy", cn2, ignore.case = TRUE)]))
  }
  stop("未知样本命名: ", paste(head(cn2), collapse = ", "))
}

score_gsva <- function(mat, gsets){
  library(GSVA)
  gsva(gsvaParam(exprData = as.matrix(mat), geneSets = gsets))
}

## ---------- 1. 基因集定义 ----------
lectin_sets <- list(
  Galectins    = c("LGALS1","LGALS2","LGALS3","LGALS4","LGALS7","LGALS8","LGALS9",
                   "LGALS12","LGALS13","LGALS14","LGALS16","LGALS3BP"),
  Siglecs      = c("SIGLEC1","SIGLEC5","SIGLEC7","SIGLEC9","SIGLEC10","SIGLEC11","SIGLEC14","SIGLEC16"),
  C_type_lectin= c("CLEC7A","CLEC4A","CLEC5A","CLEC12A","CLEC4E","MRC1","CD207"),
  Selectins    = c("SELE","SELL","SELP","SELPLG")
)

glyco_enzyme_sets <- list(
  Glycosyltransferases = c("MGAT1","MGAT2","MGAT3","MGAT4A","MGAT4B","MGAT5","MGAT5B",
                           "ST6GAL1","ST6GAL2","ST3GAL1","ST3GAL2","ST3GAL3","ST3GAL4","ST3GAL5",
                           "B4GALT1","B4GALT2","B4GALT5","B3GNT2","FUT8","FUT2"),
  Glycosidases         = c("NEU1","NEU2","NEU3","NEU4","HEXA","HEXB","GLB1","MAN2A1","MAN2B1"),
  O_GlcNAc_machine     = c("OGT","OGA","GFPT1","GFPT2","GNPNAT1","UAP1"),
  Ganglioside_synth    = c("UGCG","ST3GAL5","B4GALNT1","ST8SIA1","ST8SIA5","B3GALT4"),
  Sulfation            = c("CHST2","CHST3","CHST7","GAL3ST1","SULF1","SULF2")
)

ldam_sig  <- c("CD9","CD44","LGALS3","APOE","ABCA1","PLIN3","C1QA","C1QB","C1QC","TREM2","SPP1","LPL")
ferro_sig <- c("GPX4","SLC7A11","ACSL4","LPCAT3","ALOX15","TFRC","FTH1","FTL","CP","GCLM","GSTP1","NCOA4","HSPB1")

## ============================================================
## 补14-1. 凝集素家族全景热图（bulk + 血液 + 小鼠DEG + 蛋白组）
## ============================================================
tryCatch({
  cat("\n===== 14-1 凝集素家族全景 =====\n")
  all_lectins <- unique(unlist(lectin_sets))
  fam <- rep(names(lectin_sets), sapply(lectin_sets, length))
  names(fam) <- unlist(lectin_sets)
  fc_tab <- data.frame(gene = all_lectins)

  mats <- list(bulk_hippo = "01_GSE134697_vst_normalized.csv",
               blood      = "08_GSE143272_normalized.csv")
  for(nm in names(mats)){
    m <- read_mat(file.path(OUT_DIR, mats[[nm]]))
    g <- split_grp(colnames(m))
    cat(sprintf("  %s: 疾病 %d 例, 对照 %d 例\n", nm, length(g$dis), length(g$con)))
    hit <- intersect(all_lectins, rownames(m))
    l2 <- log2(rowMeans(m[hit, g$dis, drop=FALSE]+0.5)) - log2(rowMeans(m[hit, g$con, drop=FALSE]+0.5))
    fc_tab[[nm]] <- round(l2[match(all_lectins, hit)], 3)
  }

  ## 小鼠 DEG 表
  deg <- fread(file.path(OUT_DIR, "08_GSE88992_mouse_DEG_KA_vs_Saline.csv"), data.table = FALSE)
  cat("小鼠 DEG 列名:\n"); print(names(deg))
  gcol <- names(deg)[grep("gene|symbol", names(deg), ignore.case = TRUE)][1]
  fcol <- names(deg)[grep("log.*fc|logfc", names(deg), ignore.case = TRUE)][1]
  stopifnot(!is.na(gcol), !is.na(fcol))
  fc_tab$mouse <- round(deg[[fcol]][match(tolower(all_lectins), tolower(deg[[gcol]]))], 3)

  fwrite(fc_tab, file.path(OUT_DIR, "14_lectin_landscape_log2FC.csv"))
  print(fc_tab)

  hm <- as.matrix(fc_tab[, -1, drop = FALSE]); rownames(hm) <- fc_tab$gene
  hm <- hm[rowSums(!is.na(hm)) > 0, , drop = FALSE]
  lim <- max(abs(hm), na.rm = TRUE)
  pdf(file.path(OUT_DIR, "14_lectin_landscape_heatmap.pdf"), width = 6.5, height = 9)
  pheatmap(hm, cluster_rows = TRUE, cluster_cols = FALSE, na_col = "grey90",
           color = colorRampPalette(c("#2166AC","white","#C0392B"))(100),
           breaks = seq(-lim, lim, length.out = 101),
           annotation_row = data.frame(Family = fam[rownames(hm)], row.names = rownames(hm)),
           main = "Lectin family remodelling (log2FC, disease vs control)")
  dev.off()

  rat <- read_excel(file.path(BASE_DIR, "Differentially_expressed_statistics.xlsx"), sheet = "KAvsnc")
  raw_m <- read_excel(file.path(BASE_DIR, "12974_2023_2840_MOESM3_ESM.xlsx"), col_names = FALSE)
  hdr <- which(apply(raw_m, 1, function(x) any(grepl("^Gene name$", x)) && any(grepl("Ratio", x))))[1]
  mou <- read_excel(file.path(BASE_DIR, "12974_2023_2840_MOESM3_ESM.xlsx"), skip = hdr - 1)
  prot <- bind_rows(
    data.frame(gene = rat[["Gene name"]], FC = as.numeric(rat[["KA/nc Ratio"]]),
               P = as.numeric(rat[["KA/nc P value"]]), source = "Rat proteome"),
    data.frame(gene = mou[["Gene name"]], FC = as.numeric(mou[["TLE/Con Ratio"]]),
               P = as.numeric(mou[["TLE/Con P value"]]), source = "Mouse proteome"))
  prot <- prot[toupper(prot$gene) %in% all_lectins, ]
  fwrite(prot, file.path(OUT_DIR, "14_lectin_proteome_hits.csv"))
  print(prot)

  cat("14-1 完成\n")
}, error = function(e) cat("14-1 出错:", conditionMessage(e), "\n"))

## ============================================================
## 补14-2. LGALS3 共表达 + 桥梁定量检验
## ============================================================
tryCatch({
  cat("\n===== 14-2 LGALS3 共表达网络 =====\n")
  m <- read_mat(file.path(OUT_DIR, "01_GSE134697_vst_normalized.csv"))
  g <- split_grp(colnames(m)); dis <- g$dis
  stopifnot("LGALS3" %in% rownames(m))

  targets <- unique(c(ldam_sig, ferro_sig, "SYNGR2","LPCAT2","SPI1","NFKB1","FOS","RFX5"))
  targets <- intersect(setdiff(targets, "LGALS3"), rownames(m))
  cors <- sapply(targets, function(g2){
    ct <- suppressWarnings(cor.test(as.numeric(m["LGALS3", dis]), as.numeric(m[g2, dis]), method = "spearman"))
    c(rho = unname(ct$estimate), P = ct$p.value)
  })
  cor_df <- data.frame(gene = targets, rho = round(cors["rho",],3), P = signif(cors["P",],3),
                       axis = ifelse(targets %in% ferro_sig, "Ferroptosis", "LDAM/TF"))
  cor_df <- cor_df[order(-cor_df$rho), ]
  fwrite(cor_df, file.path(OUT_DIR, "14_LGALS3_coexpression.csv"))
  print(cor_df)

  p <- ggplot(cor_df, aes(reorder(gene, rho), rho, fill = axis)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("LDAM/TF"="#C0392B","Ferroptosis"="#2166AC")) +
    labs(x = NULL, y = "Spearman rho with LGALS3 (TLE hippocampus)",
         title = "LGALS3 co-expression bridges LDAM and ferroptosis programs") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT_DIR, "14_LGALS3_coexpression.pdf"), p, width = 7, height = 6)
  ggsave(file.path(OUT_DIR, "14_LGALS3_coexpression.png"), p, width = 7, height = 6, dpi = 300)

  sc <- score_gsva(m[, dis], list(LDAM = ldam_sig, Ferroptosis = ferro_sig))
  lg <- as.numeric(m["LGALS3", dis])
  b1 <- cor.test(lg, as.numeric(sc["LDAM",]), method = "spearman")
  b2 <- cor.test(lg, as.numeric(sc["Ferroptosis",]), method = "spearman")
  b3 <- cor.test(as.numeric(sc["LDAM",]), as.numeric(sc["Ferroptosis",]), method = "spearman")
  bridge <- data.frame(pair = c("LGALS3~LDAM","LGALS3~Ferroptosis","LDAM~Ferroptosis"),
                       rho = round(c(b1$estimate, b2$estimate, b3$estimate), 3),
                       P = signif(c(b1$p.value, b2$p.value, b3$p.value), 3))
  print(bridge); fwrite(bridge, file.path(OUT_DIR, "14_LGALS3_bridge_test.csv"))
  cat("14-2 完成（重点看 bridge 表三个 rho）\n")
}, error = function(e) cat("14-2 出错:", conditionMessage(e), "\n"))

## ============================================================
## 补14-3. 糖基化酶功能分组打分（bulk + 单细胞）完整版
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr)
})

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")

## ---------- 工具 ----------
read_mat <- function(path){
  d <- fread(path, data.table = FALSE)
  m <- as.matrix(d[, -1, drop = FALSE])
  rownames(m) <- make.unique(as.character(d[[1]]))
  m
}
split_grp <- function(cn){
  cn2 <- gsub("\\.", "-", cn)
  if(all(grepl("^CHKJ", cn2))){
    return(list(dis = cn[grepl("^CHKJ-H", cn2)], con = cn[grepl("^CHKJ-C", cn2)]))
  }
  if(any(grepl("Healthy", cn2))){
    return(list(dis = cn[grepl("Drug", cn2, ignore.case = TRUE)],
                con = cn[grepl("Healthy", cn2, ignore.case = TRUE)]))
  }
  stop("未知样本命名: ", paste(head(cn2), collapse = ", "))
}
score_gsva <- function(mat, gsets){
  library(GSVA)
  gsva(gsvaParam(exprData = as.matrix(mat), geneSets = gsets))
}

## ---------- 基因集 ----------
glyco_enzyme_sets <- list(
  Glycosyltransferases = c("MGAT1","MGAT2","MGAT3","MGAT4A","MGAT4B","MGAT5","MGAT5B",
                           "ST6GAL1","ST6GAL2","ST3GAL1","ST3GAL2","ST3GAL3","ST3GAL4","ST3GAL5",
                           "B4GALT1","B4GALT2","B4GALT5","B3GNT2","FUT8","FUT2"),
  Glycosidases         = c("NEU1","NEU2","NEU3","NEU4","HEXA","HEXB","GLB1","MAN2A1","MAN2B1"),
  O_GlcNAc_machine     = c("OGT","OGA","GFPT1","GFPT2","GNPNAT1","UAP1"),
  Ganglioside_synth    = c("UGCG","ST3GAL5","B4GALNT1","ST8SIA1","ST8SIA5","B3GALT4"),
  Sulfation            = c("CHST2","CHST3","CHST7","GAL3ST1","SULF1","SULF2")
)

## ============================================================
## A. bulk 海马（GSE134697）
## ============================================================
tryCatch({
  cat("\n===== 14-3A bulk 糖基化酶谱 =====\n")
  m <- read_mat(file.path(OUT_DIR, "01_GSE134697_vst_normalized.csv"))
  g <- split_grp(colnames(m))
  cat(sprintf("bulk: 疾病 %d 例, 对照 %d 例\n", length(g$dis), length(g$con)))
  sc <- score_gsva(m, glyco_enzyme_sets)
  st <- t(sapply(rownames(sc), function(gs){
    w <- wilcox.test(sc[gs, g$dis], sc[gs, g$con])
    c(delta = mean(sc[gs, g$dis]) - mean(sc[gs, g$con]), P = w$p.value)
  }))
  enz <- data.frame(program = rownames(st), delta = round(st[,"delta"],3), P = signif(st[,"P"],3))
  print(enz); fwrite(enz, file.path(OUT_DIR, "14_glyco_enzyme_bulk.csv"))

  p <- ggplot(enz, aes(reorder(program, delta), delta, fill = P < 0.05)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("TRUE"="#C0392B","FALSE"="grey70")) +
    labs(x = NULL, y = "Score difference (TLE - control)",
         title = "Glycosylation enzyme programs in epileptic hippocampus", fill = "P < 0.05") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT_DIR, "14_glyco_enzyme_bulk.pdf"), p, width = 7, height = 4.5)
  ggsave(file.path(OUT_DIR, "14_glyco_enzyme_bulk.png"), p, width = 7, height = 4.5, dpi = 300)
  cat("14-3A 完成\n")
}, error = function(e) cat("14-3A 出错:", conditionMessage(e), "\n"))

## ============================================================
## 补14-3B + 14-4 完整修正版 v2
## ============================================================
suppressPackageStartupMessages({ library(CellChat); library(UCell); library(data.table); library(dplyr) })

if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root
OUT_DIR  <- file.path(BASE_DIR, "results")

glyco_enzyme_sets <- list(
  Glycosyltransferases = c("MGAT1","MGAT2","MGAT3","MGAT4A","MGAT4B","MGAT5","MGAT5B",
                           "ST6GAL1","ST6GAL2","ST3GAL1","ST3GAL2","ST3GAL3","ST3GAL4","ST3GAL5",
                           "B4GALT1","B4GALT2","B4GALT5","B3GNT2","FUT8","FUT2"),
  Glycosidases         = c("NEU1","NEU2","NEU3","NEU4","HEXA","HEXB","GLB1","MAN2A1","MAN2B1"),
  O_GlcNAc_machine     = c("OGT","OGA","GFPT1","GFPT2","GNPNAT1","UAP1"),
  Ganglioside_synth    = c("UGCG","ST3GAL5","B4GALNT1","ST8SIA1","ST8SIA5","B3GALT4"),
  Sulfation            = c("CHST2","CHST3","CHST7","GAL3ST1","SULF1","SULF2")
)

cc0 <- readRDS(file.path(OUT_DIR, "06_CellChat_object.rds"))
expr <- cc0@data
mg   <- cc0@meta

## ============================================================
## 14-3B 单细胞糖基化酶谱
## ============================================================
tryCatch({
  cat("\n===== 14-3B 单细胞糖基化酶谱 =====\n")
  keep_genes <- intersect(unique(unlist(glyco_enzyme_sets)), rownames(expr))
  cat(sprintf("打分基因: %d/%d 检出\n", length(keep_genes), length(unique(unlist(glyco_enzyme_sets)))))

  sc_raw <- ScoreSignatures_UCell(expr[keep_genes, ], features = glyco_enzyme_sets)
  ## 自适应方向：哪一维等于细胞数，哪一维就是细胞
  if(nrow(sc_raw) == ncol(expr)){
    scores <- as.matrix(sc_raw)          # 已是 细胞 x 程序
  } else if(ncol(sc_raw) == ncol(expr)){
    scores <- t(as.matrix(sc_raw))       # 程序 x 细胞，转置
  } else stop("UCell 输出维度无法识别: ", paste(dim(sc_raw), collapse=" x "))
  rownames(scores) <- colnames(expr)
  cat("scores:", nrow(scores), "细胞 x", ncol(scores), "程序\n")

  v0 <- mg[["LDAM_level"]]
  hi <- grepl("high", as.character(v0), ignore.case = TRUE)
  hi[is.na(v0)] <- NA
  cat("LDAM-high:", sum(hi, na.rm = TRUE), " / LDAM-low:", sum(!hi, na.rm = TRUE), "\n")
  stopifnot(identical(rownames(mg), colnames(expr)))

  out <- lapply(colnames(scores), function(pr){
    v <- scores[, pr]
    keep <- !is.na(hi) & !is.na(v)
    w <- wilcox.test(v[keep & hi], v[keep & !hi])
    data.frame(program = pr,
               delta = round(mean(v[keep & hi]) - mean(v[keep & !hi]), 4),
               P = signif(w$p.value, 3))
  })
  sc_enz <- do.call(rbind, out)
  print(sc_enz)
  fwrite(sc_enz, file.path(OUT_DIR, "14_glyco_enzyme_singlecell.csv"))

  if("LGALS3" %in% rownames(expr)){
    lg3 <- as.numeric(expr["LGALS3", ])
    keep <- !is.na(hi)
    w <- wilcox.test(lg3[keep & hi], lg3[keep & !hi])
    cat(sprintf("\n[关键] 单细胞 LGALS3: LDAM-high mean=%.4f vs LDAM-low mean=%.4f, 检出率 %.1f%% vs %.1f%%, P=%.3e\n",
                mean(lg3[keep & hi]), mean(lg3[keep & !hi]),
                100*mean(lg3[keep & hi]>0), 100*mean(lg3[keep & !hi]>0), w$p.value))
  }
  cat("14-3B 完成\n")
}, error = function(e) cat("14-3B 出错:", conditionMessage(e), "\n"))

## ============================================================
## 14-4 LGALS3 细胞通讯（chat_group 含 LDAM-high/low 小胶质细胞）
## ============================================================
tryCatch({
  cat("\n===== 14-4 LGALS3 细胞通讯 =====\n")
  cat("chat_group 水平:", paste(unique(mg$chat_group), collapse = ", "), "\n")

  gal_lr <- data.frame(
    interaction_name = paste0("LGALS3_", c("TREM2","CD44","ITGB1","BSG","LRP1","SLC3A2","PTPRC")),
    pathway_name     = "GALECTIN",
    ligand   = "LGALS3",
    receptor = c("TREM2","CD44","ITGB1","BSG","LRP1","SLC3A2","PTPRC"),
    agonist = "", antagonist = "", co_A_receptor = "", co_I_receptor = "",
    evidence = "Literature-curated Galectin-3 binding partners",
    annotation = "Secreted Signaling",
    interaction_name_2 = paste("LGALS3 -", c("TREM2","CD44","ITGB1","BSG","LRP1","SLC3A2","PTPRC")),
    stringsAsFactors = FALSE)
  db <- CellChatDB.human
  db$interaction <- bind_rows(db$interaction, gal_lr)   # bind_rows 自动对齐/补齐列

  cc <- createCellChat(object = expr, meta = mg, group.by = "chat_group")
  cc@DB <- db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cat("计算通讯概率（62k 细胞，耐心等几分钟）...\n")
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10)

  lr <- subsetCommunication(cc)
  cat("全部互作对:", nrow(lr), "\n")
  lg <- lr[grepl("LGALS", lr$ligand) | grepl("LGALS", lr$receptor), ]
  cat("含 LGALS 的互作对:", nrow(lg), "\n")
  if(nrow(lg) > 0){
    fwrite(lg, file.path(OUT_DIR, "14_CellChat_LGALS3_pairs.csv"))
    show_cols <- intersect(c("source","target","ligand","receptor","prob","pval"), names(lg))
    lg <- lg[order(-lg$prob), ]
    print(lg[, show_cols])

    ## 核心问题：LDAM-high vs LDAM-low 小胶质细胞作为信号"发出者"的对比
    send <- lg[lg$ligand == "LGALS3", ]
    tab <- aggregate(prob ~ source + receptor, send, sum)
    cat("\n--- LGALS3 发出强度按来源细胞汇总 ---\n"); print(tab[order(-tab$prob), ])

    pdf(file.path(OUT_DIR, "14_CellChat_LGALS3_bubble.pdf"), width = 10, height = 5.5)
    print(netVisual_bubble(cc, pairLR.use = unique(lg[, c("interaction_name","pathway_name")]),
                           remove.isolate = TRUE, title.name = "LGALS3-mediated signalling"))
    dev.off()
  } else {
    cat("仍为 0——检查 LGALS3 与受体基因在数据中的检出率：\n")
    for(g in c("LGALS3","TREM2","CD44","ITGB1","BSG","LRP1","SLC3A2","PTPRC")){
      if(g %in% rownames(expr)) cat(sprintf("  %-8s 检出率 %.2f%%\n", g, 100*mean(expr[g, ]>0)))
      else cat(sprintf("  %-8s 不在矩阵\n", g))
    }
  }
  saveRDS(cc, file.path(OUT_DIR, "14_CellChat_gal_custom.rds"))
  cat("14-4 完成（对象已存盘）\n")
}, error = function(e) cat("14-4 出错:", conditionMessage(e), "\n"))

cat("\n全部完成\n")
## ============================================================
## 补14-5. 血液 LGALS3 + 凝集素 signature（GSE143272）
## ============================================================
tryCatch({
  cat("\n===== 14-5 血液 LGALS3 =====\n")
  m <- read_mat(file.path(OUT_DIR, "08_GSE143272_normalized.csv"))
  g <- split_grp(colnames(m))
  cat(sprintf("血液: 疾病 %d 例, 对照 %d 例\n", length(g$dis), length(g$con)))
  stopifnot("LGALS3" %in% rownames(m))

  w <- wilcox.test(as.numeric(m["LGALS3", g$dis]), as.numeric(m["LGALS3", g$con]))
  cat(sprintf("血液 LGALS3: 疾病 mean=%.3f, 对照 mean=%.3f, Wilcoxon P=%.4f\n",
              mean(m["LGALS3", g$dis]), mean(m["LGALS3", g$con]), w$p.value))

  sc <- score_gsva(m, list(Galectins = lectin_sets$Galectins,
                           Lectins_all = unique(unlist(lectin_sets))))
  bl <- data.frame(sample = colnames(m),
                   group = ifelse(colnames(m) %in% g$dis, "TLE", "Control"),
                   t(sc), check.names = FALSE)
  bl_long <- pivot_longer(bl, -c(sample, group))
  p <- ggplot(bl_long, aes(group, value, fill = group)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, alpha = .5) +
    facet_wrap(~name, scales = "free") +
    scale_fill_manual(values = c("TLE"="#C0392B","Control"="#2166AC")) +
    labs(title = "Blood lectin signatures (GSE143272)", y = "ssGSEA score") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(OUT_DIR, "14_blood_lectin.pdf"), p, width = 8, height = 4.5)
  fwrite(bl, file.path(OUT_DIR, "14_blood_lectin_scores.csv"))
  cat("14-5 完成\n")
}, error = function(e) cat("14-5 出错:", conditionMessage(e), "\n"))

cat("\n全部完成，结果在: ", OUT_DIR, "\n")
