# ============================================================================
# Script      : 跨组学蛋白验证.R
# Purpose     : Cross-omics protein-level validation: joint reading of rat KA and mouse chronic TLE proteome tables with header auto-detection.
# Manuscript  : Supplementary analyses (cross-species proteome)
# Input       : Differentially_expressed_statistics.xlsx; 12974_2023_2840_MOESM3_ESM.xlsx
# Output      : cross-omics comparison tables/figures
# Dependencies: dplyr, ggplot2, readxl
# Date        : standardized 2026-10 (original analyses 2025-09/10)
# Paths       : absolute paths replaced by BASE_DIR = Sys.getenv("LDAM_BASE", ".")
#               (see _config.R; run from the project root or set LDAM_BASE)
# ============================================================================

## ---- 1. 读入两个蛋白组（自动定位表头） ----
if (!exists("BASE_DIR")) BASE_DIR <- Sys.getenv("LDAM_BASE", ".")  # project root: set env LDAM_BASE or run from project root

library(readxl); library(dplyr); library(ggplot2)

smart_read <- function(path, sheet = 1){
  raw <- read_excel(path, sheet = sheet, col_names = FALSE)
  hdr <- which(apply(raw, 1, function(x)
    any(grepl("^Gene name$", x, ignore.case = TRUE)) &&
      any(grepl("Ratio", x, ignore.case = TRUE))))[1]
  if (is.na(hdr)) stop("找不到真表头行: ", path)
  df <- read_excel(path, sheet = sheet, skip = hdr - 1)
  df <- df[, !is.na(names(df)) & names(df) != "" & !grepl("^\\.\\.\\.", names(df))]
  df
}

rat   <- smart_read(file.path(BASE_DIR, "Differentially_expressed_statistics.xlsx"), sheet = "KAvsnc")
mouse <- smart_read(file.path(BASE_DIR, "12974_2023_2840_MOESM3_ESM.xlsx"), sheet = 1)

cat("大鼠关键列:\n");  print(grep("Gene|Ratio|P value", names(rat), value = TRUE))
cat("小鼠关键列:\n");  print(grep("Gene|Ratio|P value", names(mouse), value = TRUE))

## 按实际列名提取（模糊匹配，防止斜杠/大小写差异）
grab <- function(df, key) names(df)[grepl(key, names(df), ignore.case = TRUE)][1]

df <- bind_rows(
  data.frame(gene = rat[[grab(rat, "Gene name")]],
             FC   = as.numeric(rat[[grab(rat, "KA/nc Ratio")]]),
             P    = as.numeric(rat[[grab(rat, "KA/nc P value")]]),
             source = "Rat KA proteome (this study)"),
  data.frame(gene = mouse[[grab(mouse, "Gene name")]],
             FC   = as.numeric(mouse[[grab(mouse, "TLE/Con Ratio")]]),
             P    = as.numeric(mouse[[grab(mouse, "TLE/Con P value")]]),
             source = "Mouse TLE proteome (published)")
)

## ---- 2. 目标蛋白 ----
targets <- list(
  "LDAM signature"   = c("C1qa","C1qb","C1qc","Cd44","Lgals3","Cd9","Apoe","Abca1","Plin3"),
  "Hub genes"        = c("Syngr2","Hspb1"),
  "Ferroptosis/iron" = c("Cp","Ftl1","Gclm","Gstp1","Por","Glud1")
)
gene2cat <- setNames(rep(names(targets), sapply(targets, length)), tolower(unlist(targets)))

df <- df %>%
  filter(!is.na(FC), !is.na(P), P < 0.05, FC > 1) %>%
  mutate(key = tolower(gene)) %>%
  filter(key %in% names(gene2cat)) %>%
  mutate(gene = toupper(gene), category = gene2cat[key])
df$key <- NULL   # 代替 select(-key)

cat("命中记录数:", nrow(df), "\n"); print(df)
write.csv(df, file.path(OUT_DIR, "13_crossomics_proteome_hits.csv"), row.names = FALSE)

## ---- 3. 森林图 ----
p <- ggplot(df, aes(x = log2(FC), y = reorder(gene, FC), color = source, shape = source)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  scale_color_manual(values = c("Rat KA proteome (this study)" = "#C0392B",
                                "Mouse TLE proteome (published)" = "#2166AC")) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(x = "log2 fold change (protein level)", y = NULL, color = NULL, shape = NULL,
       title = "Cross-omics validation of the LDAM-ferroptosis axis") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.text.y = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "13_crossomics_forest.pdf"), p, width = 7.5, height = 8)
ggsave(file.path(OUT_DIR, "13_crossomics_forest.png"), p, width = 7.5, height = 8, dpi = 300)
cat("13_跨组学蛋白验证 完成\n")
