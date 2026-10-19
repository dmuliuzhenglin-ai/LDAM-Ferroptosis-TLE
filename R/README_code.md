# 代码清单与运行说明（README_code）

本目录为手稿配套 R 分析脚本的**投稿级标准化副本**。原脚本（项目根目录 `R/` 及根目录下 `10.R`、`11_高级补充分析.R`）未做任何改动。

## 标准化内容

- 每个脚本顶部加入统一英文头注释块（Script / Purpose / Manuscript / Input / Output / Dependencies / Date: 2026-10）。
- 全部硬编码绝对路径（`C:/Users/...` 共 58 处）替换为项目根变量：`BASE_DIR <- Sys.getenv("LDAM_BASE", ".")`，统一定义见 `_config.R`；每个脚本自带守护行，从项目根目录运行或设置环境变量 `LDAM_BASE` 即可。
- 中文注释保留；主要中文段落标题下补充一行英文摘要（`# [EN] ...`）。
- 调试残留处理：纯查看型 `print(names()/dim()/head()...)` 已注释并标 `# [debug]`；一次性装包语句（`install.packages` / `BiocManager::install` / `remotes::install_github`）已注释并标 `# [setup]`（首次运行请先跑 `00_setup_genesets.R`）。
- 编码统一 UTF-8（LF 换行），行尾空白已清除，原缩进风格保留（安全第一，未做重排）。

## 敏感信息处理

- `15_MR分析.R` 第 1 行原含一枚 **OpenGWAS JWT 访问令牌**（RS256，payload 内含账号邮箱），已删除并替换为环境变量读取 + 显式报错提示。**建议令牌持有者尽快到 OpenGWAS 后台吊销/轮换该令牌**。
- 全部脚本中的本机用户名路径（`C:/Users/19367/...`）已随路径参数化一并移除。
- 其余 `secret`/`token` 关键词匹配均为基因集名（如 `Rab_secretory`）等误报，无其他敏感信息。

## 遗留与异常文件说明

- `07_scRNA_GSE201048_validation.R`、`fix_04前置修复.R`：内容与文件名不符，为早期胶质瘤项目遗留代码，与本手稿无关；已标注 LEGACY/DEPRECATED，外部路径参数化为 `GLIOMA_GEO_DIR` / `GLIOMA_TCGA_DIR` 环境变量。
- `微胶质细胞 pseudobulk 偏相关.R`：源文件含一段旧版函数体的孤儿残留（原文件本身无法通过 R 解析），已整段注释并标注 `# [deprecated]`，不影响主流程。
- `37_blood_combined_ROC_v2.R` 取代 `37_blood_combined_ROC.R`；`28修正.R` 取代 `28_sc_exosome_modules.R`（旧版均保留备查）。
- 空占位文件 `R/23外泌体.R`（0 字节）未收录；项目根目录 `10.R` 与 `R/10.R` 逐字节相同，只收录一份。

## R 版本声明

脚本 developed under R 4.x（原脚本使用 R ≥ 4.0 的 raw string 语法及 Seurat v5 API，建议 R ≥ 4.3 + 最新 Bioconductor）。本机无 Rscript，**未做 `Rscript -e 'parse(...)'` 运行级验证**；替代验证为文本级括号/引号配平检查（字符串感知，含 R raw string 与多行字符串），全部 46 个文件通过。

## 建议运行顺序

```
00_setup_genesets.R            # 第 1 个运行，只需一次（装包 + 冻结基因集）
01 -> 02 -> 02b -> 03          # bulk 发现 + 验证 + WGCNA
08 (Part A)                    # 外周血预处理，04 依赖其输出
04 -> 05                       # 机器学习 hub 基因 + 免疫/富集
06                             # 核心单细胞（后续多数单细胞脚本依赖其 rds）
补充.R -> 28修正.R -> 30 -> 31 -> 34 -> 36 -> 微胶质细胞 pseudobulk 偏相关.R
15 -> 16 -> 17 -> 38 -> FigureS5.R        # 遗传证据链（MR/SMR/coloc）
扰动拯救分析.R -> 35.R                    # 扰动/拯救
Gal-3 受体轴分析.R -> NicheNet 配体-靶标分析.R
糖基化.R -> 14_糖基化补充分析.R
23 -> 27 -> 27b -> 40 -> 跨组学蛋白验证.R  # 外泌体/蛋白组
24 -> 25 -> 26 -> 29                      # 转录后/miRNA
32                                          # 药物重定位
33 -> 37_blood_combined_ROC_v2.R            # 血液 ROC
10.R                                        # 最后出主图
11_高级补充分析.R                            # 探索性补充，任意时间
```

## 脚本清单

| 文件 | 分析内容 | 手稿章节 | 依赖 R 包 |
|---|---|---|---|
| `00_setup_genesets.R` | 环境自检/装包 + 构建并冻结 LDAM 与铁死亡基因集（只跑一次） | Methods 2.1 (Datasets and gene sets) | data.table, dplyr, ggplot2, msigdbr, stringr |
| `01_bulk_GSE134697_paired_DE.R` | GSE134697 配对 DESeq2 差异分析（17 例 mTLE 海马 vs 自体新皮层） | Methods 2.2 / Results 3.1 | DESeq2, data.table, dplyr, ggplot2, ggrepel, org.Hs.eg.db, stringr |
| `02_bulk_validation_cohorts_DE.R` | 验证队列差异分析：GSE190451（人 TLE 皮层 bulk）+ GSE168375（纯化星形胶质细胞） | Methods 2.2 / Results 3.1 | DESeq2, GEOquery, data.table, dplyr, ggplot2, ggrepel, limma, org.Hs.eg.db, stringr |
| `02b_GSE190451_fix.R` | 补跑 GSE190451 并与 GSE168375 做 log2FC 方向一致性检查 | Methods 2.2 / Results 3.1 | DESeq2, data.table, ggplot2, limma |
| `03_bulk_GSE63808_WGCNA.R` | GSE63808（129 例 TLE 海马）WGCNA 共表达模块 + LDAM/铁死亡 ssGSEA 关联 | Methods 2.2 / Results 3.3 | GEOquery, GSVA, WGCNA, data.table, ggplot2, limma, stringr |
| `04_machine_learning_hub_ROC.R` | LASSO + SVM-RFE + 随机森林三算法交集筛 hub 基因；组织与外周血 ROC | Methods 2.3 / Results 3.2 | caret, data.table, dplyr, e1071, ggplot2, glmnet, pROC, randomForest |
| `05_immune_GSVA_enrichment.R` | GO/KEGG/GSEA 富集 + ssGSEA 免疫浸润 + GSVA LDAM/铁死亡评分相关 | Methods 2.2 / Results 3.3 | GSVA, clusterProfiler, data.table, dplyr, enrichplot, ggplot2, ggpubr, org.Hs.eg.db, pheatmap |
| `06_scRNA_GSE190452_main.R` | 核心 snRNA-seq（GSE190452，4v4）：整合、注释、LDAM/铁死亡打分、CellChat | Methods 2.4 / Results 3.4-3.6 | AUCell, CellChat, Seurat, SeuratWrappers, UCell, data.table, dplyr, ggplot2, ggpubr, harmony, monocle3, patchwork, stringr |
| `07_scRNA_GSE201048_validation.R` | 【遗留模板】旧胶质瘤项目 BayesPrism 去卷积流程，与本手稿无关，仅供留档 | not used in manuscript (legacy file, kept for transparency) | BayesPrism, Matrix, Seurat, data.table, survival |
| `08_blood_GSE143272_mouse_GSE88992.R` | GSE143272 外周血预处理（供 04 做血液 ROC）+ GSE88992 小鼠 KA 模型旁证 | Methods 2.5 / Results 3.7 | GSVA, affy, biomaRt, limma, mouse4302.db |
| `10.R` | 出版级主图拼装：从 results/ 直接出 7 张主图（PDF+PNG） | Main Figures 1-7 | CellChat, Seurat, data.table, dplyr, ggplot2, pROC, patchwork, tidyr |
| `11_高级补充分析.R` | 高级补充探索：拟时序校验、调控子、共识聚类、miRNA 关联（探索性） | Supplementary analyses | CellChat, ConsensusClusterPlus, GSVA, Seurat, SeuratWrappers, data.table, decoupleR, dplyr, ggplot2, ggraph, igraph, monocle3, multiMiR, pROC, tidyr |
| `14_糖基化补充分析.R` | 糖基化补充：凝集素全景、LGALS3 共表达/桥梁、糖基化酶谱、CellChat、血液 | Supplementary analyses (Galectin-3 axis related) | CellChat, GSVA, UCell, data.table, dplyr, ggplot2, pheatmap, readxl, tidyr |
| `15_MR分析.R` | 双样本 MR：Gal-3 pQTL(deCODE)/血脂/hub 基因 eQTL → 癫痫 GWAS（OpenGWAS） | Methods 2.10 / Results 3.11 | TwoSampleMR, coloc, data.table, dplyr, ieugwasr |
| `16_脑组织eQTL_MR.R` | LGALS3 脑组织 eQTL（MetaBrain 海马/皮层）→ 癫痫结局 MR | Methods 2.10 / Results 3.11 | TwoSampleMR, data.table, dplyr, ieugwasr |
| `17_SMR_HEIDI敏感性分析.R` | LGALS3 皮层 eQTL × 局灶伴 HS 癫痫 GWAS 的 SMR + HEIDI 敏感性分析 | Methods 2.12 / Results 3.13 | data.table, dplyr |
| `23_exosome_proteome_enrichment.R` | 大鼠 KA 海马 TMT 蛋白组中外泌体蛋白富集 + Gal-3/CD9 囊泡定位 + 跨物种重复 | Supplementary analyses (cross-species proteome) | dplyr, ggplot2, readr, stringr, tibble, tidyr |
| `24.R` | 转录后闸门检验：LGALS3 蛋白-mRNA 不一致、mTORC1/翻译、蛋白酶体/自噬、外泌体装载 | Supplementary analyses (post-transcriptional regulation) | dplyr, ggplot2, readr, stringr, tibble, tidyr |
| `25_miRNA_LGALS3_TLE.R` | 靶向 LGALS3 的 miRNA 在 TLE 海马是否下调（GSE99455 主队列 + GSE245228 验证） | Supplementary analyses (miRNA axis) | dplyr, ggplot2, readr, stringr, tibble |
| `26_miR34a_target_derepression.R` | miR-34a-5p 靶基因集体去抑制检验（通路级证据） | Supplementary analyses (miRNA axis) | dplyr, ggplot2, readr, stringr, tibble |
| `27_exosome_function.R` | 外泌体机器（摄取/装载/分泌）在人慢性期与小鼠 KA 时序中的功能分析 | Supplementary analyses (extracellular-vesicle modules) | ggplot2, tidyverse |
| `27b_exosome_mouse_proteome.R` | 外部独立验证：慢性小鼠 TLE 海马蛋白组（KA 35 天 LFQ）外泌体模块富集 | Supplementary analyses (cross-species proteome) | readxl, tidyverse |
| `28_sc_exosome_modules.R` | 单细胞外泌体模块 UCell 评分：LDAM-high vs LDAM-low 小胶质细胞 | Supplementary analyses (single-cell EV modules) | Seurat, UCell, ggplot2, tidyverse |
| `28修正.R` | 28 的 v2 修正版（UCell 自动后缀处理），取代 28_sc_exosome_modules.R | Supplementary analyses (single-cell EV modules) | Seurat, UCell, ggplot2, tidyverse |
| `29_miRNA_derepression_panel.R` | 其余候选 miRNA 的靶基因集体去抑制面板检验（TargetScan 8.0） | Supplementary analyses (miRNA axis) | tidyverse |
| `30_TF_activity_LDAM.R` | LDAM-high 小胶质细胞转录因子活性推断（decoupleR + DoRothEA ULM） | Supplementary analyses (upstream regulators) | Seurat, UCell, decoupleR, dorothea, ggplot2, tidyverse |
| `31_pseudotime_LDAM.R` | Monocle3 拟时序：LDAM-low → LDAM-high 轨迹及铁死亡/外泌体评分动态 | Supplementary analyses (trajectory) | Seurat, UCell, ggplot2, monocle3, tidyverse |
| `32_LINCS_drug_repositioning.R` | 基于 LDAM-high 差异签名的 LINCS L1000CDS2 药物重定位 | Methods 2.8 / Results 3.10 | Seurat, UCell, httr, jsonlite, tidyverse |
| `33_blood_exosome_scores.R` | 外泌体机器模块在外周血的可检测性 + 逐基因血液标志物候选 | Supplementary analyses (peripheral blood) | ggplot2, tidyverse |
| `34_TF_LGALS3_regulation.R` | TF 活性 → LGALS3 表达的偏相关（控制 LDAM 评分；PPARG/NR1H2/NR1H3） | Supplementary analyses (upstream regulators) | OmnipathR, Seurat, decoupleR, dplyr, ggplot2, ppcor, tibble, tidyr |
| `35.R` | Gal-3 与 DAM 程序反馈关系（JCI 2023 重分析：Tau22 / Lgals3-KO / TD139） | Supplementary analyses (perturbation, Methods 2.7 related) | dplyr, readxl, tidyr |
| `36_LGALS3_exosome_coupling.R` | LGALS3 与外泌体模块特异性共变检验（控制 LDAM，含内对照） | Supplementary analyses | Seurat, dplyr, ggplot2, ppcor, tibble, tidyr |
| `37_blood_combined_ROC.R` | 血液双 signature 联合建模 ROC（表观 + 重复 10 折 CV，DeLong） | Supplementary analyses (blood ROC; Results 3.2 related) | ggplot2, pROC, tidyverse |
| `37_blood_combined_ROC_v2.R` | 37 的完整修复版（分层分折等），取代 37_blood_combined_ROC.R | Supplementary analyses (blood ROC; Results 3.2 related) | ggplot2, pROC, tidyverse |
| `38_LGALS3_coloc_SMR.R` | LGALS3 与癫痫 GWAS 共定位（coloc.abf）+ SMR/HEIDI | Methods 2.10 / Results 3.11 | coloc, data.table, ieugwasr, tidyverse |
| `40_exosome_cargo_verification.R` | QuickGO GO:0070062 外泌体货物数据库验证（Gal-3 是否为货物蛋白等） | Supplementary analyses (EV cargo) | httr, jsonlite, tidyverse |
| `FigureS5.R` | Galectin-3 遗传证据森林图（pQTL + eQTL MR + SMR/HEIDI），Figure S5 | Figure S5 | dplyr, ggplot2, patchwork |
| `Gal-3 受体轴分析.R` | Gal-3 凝集素-糖链-受体轴四层证据（表达/共表达/CellChat/扰动），Figure S9 | Methods 2.6 / Results 3.8 (Figure S9) | Seurat, dplyr, ggplot2, ggrepel, patchwork, readxl |
| `NicheNet 配体-靶标分析.R` | NicheNet：小胶质 LGALS3(配体)→神经元靶基因预测及铁死亡富集，Figure S8 | Methods 2.6 related supplementary (Figure S8) | Seurat, dplyr, ggplot2, ggrepel, nichenetr, patchwork |
| `fix_04前置修复.R` | 【遗留】旧项目会话续跑助手（Harmony→聚类→注释），依赖会话中已有对象 | not used in manuscript (legacy helper) | Matrix, Seurat, SeuratObject, ggplot2 |
| `微胶质细胞 pseudobulk 偏相关.R` | 样本级 pseudobulk 偏相关：LGALS3×铁死亡（控制 LDAM）及反向对照，Figure S7 | Methods 2.6 / Results 3.8 (Figure S7) | Seurat, dplyr, ggplot2, ppcor |
| `扰动拯救分析.R` | 扰动拯救重分析（JCI 2023）：Lgals3 KO（Tau22 小鼠）与 TD139（人 iMGL），Figure S6 | Methods 2.7 / Results 3.9 (Figure S6) | dplyr, ggplot2, patchwork, readxl |
| `糖基化.R` | 糖基化-LDAM-铁死亡轴验证性分析（复用现有结果，离线） | Supplementary analyses (glycosylation) | AnnotationDbi, GSVA, Seurat, UCell, data.table, dplyr, ggplot2, org.Hs.eg.db, tidyr |
| `补充.R` | 单细胞补充：LDAM 高低差异基因+GO、hub 细胞类型特异性、pseudobulk 差异、小鼠 KA hub 表达 | Methods 2.4 supplementary | DESeq2, Seurat, clusterProfiler, data.table, dplyr, ggplot2, org.Hs.eg.db, pROC |
| `跨组学蛋白验证.R` | 跨组学蛋白验证：联合读取大鼠 KA 与小鼠慢性 TLE 蛋白组表（表头自动定位） | Supplementary analyses (cross-species proteome) | dplyr, ggplot2, readxl |
