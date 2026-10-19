# -*- coding: utf-8 -*-
"""
p25_deconvolution.py
====================
GSE63808 (129 TLE hippocampus bulk microarray, GPL6947) 细胞类型反卷积
参考: GSE190452 TLE 人海马 snRNA-seq (results/43_scDRS/gse190452_annotated.h5ad)

方法学要点（如实记录）:
1. bulk 表达矩阵: 解析 GSE63808_series_matrix.txt.gz (!series_matrix_table_begin 之后的
   探针×样本矩阵)；若 max>100 则 log2(x+1)（与 R/03_bulk_GSE63808_WGCNA.R 一致）。
2. 探针→基因: GPL6947.annot.gz (NCBI GEO 平台注释, 下载于分析当日)，多 symbol 取第一个，
   同基因多探针取均值（对应 limma::avereps）。
3. 签名矩阵: 单细胞 counts 层经 sc.pp.normalize_total(target_sum=1e4)+log1p 后，
   按 celltype 取平均 log-normalized 表达。
   注: 参考中 Pericyte 无注释细胞，实际使用 7 类细胞。
4. 特征基因: scanpy tl.rank_genes_groups (wilcoxon)，每类 top 100，合并去重后与
   bulk 基因级矩阵取交集。为控制 wilcoxon 运行时间，差异检验在每类 ≤3000 细胞的
   分层子抽样上进行（随机种子固定）；签名均值仍用全部细胞计算。
5. 反卷积: scipy.optimize.nnls 逐样本。先试了**绝对水平对数空间 NNLS + 截断**（bulk log2
   与签名 log1p 截断负值后直接 NNLS），结果退化（ExcNeuron 全为 0、Endothelial≈47%，
   与海马细胞构成不符，ExcNeuron 与 SNAP25 无相关），原因是 bulk 芯片与 sc 参考
   基线/量纲不可比，NNLS 主要拟合基线偏移。故主分析改用**基因中心差分对数空间 NNLS**：
   对每个特征基因，签名在 7 类细胞间减去该基因均值，bulk 在 129 样本间减去该基因均值，
   使 NNLS 拟合"相对差异"信号（NNLS 仅约束系数非负，A/b 允许负值），系数归一化
   (行和=1)为比例。敏感性分析: 线性空间 NNLS (bulk 2^x, 签名 expm1——
   expm1(均值log)≠均值exp，属近似，仅作稳健性对照)。
6. 验证: (a) 微胶质比例 vs 样本 LDAM / 铁死亡 ssGSEA 分数 Spearman;
   (b) 微胶质比例 vs bulk AIF1/CX3CR1/P2RY12 表达 Spearman (内验证);
   (c) 兴奋性神经元比例 vs SNAP25 表达 (反向 sanity check)。
7. 分层: 129 例均为 epilepsy，无对照（!Sample_characteristics 逐行核实），不做组间/预后分析。
   按估计微胶质比例 tertile 分高/低组（中 tertile 剔除），Mann-Whitney 比较两组
   bulk LDAM / 铁死亡 ssGSEA 分数；另附中位数二分组敏感性结果。

运行: gsmap_env/Scripts/python.exe scripts/p25_deconvolution.py
输出: results/63_deconv/{deconv_proportions.csv, deconv_validation.csv,
      stratification_results.csv, deconv_summary.json, signature_matrix.csv}
      figures/v25/{Fig_deconv_validation.png, Fig_deconv_stratification.png}
"""
import sys, json, gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import nnls
from scipy.stats import spearmanr, mannwhitneyu
import anndata as ad
import scanpy as sc

BASE = Path(__file__).resolve().parents[1]
RES = BASE / "results" / "63_deconv"
FIG = BASE / "figures" / "v25"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

F_SERIES = BASE / "GSE63808" / "GSE63808_series_matrix.txt.gz"
F_GPL = BASE / "GSE63808" / "GPL6947.annot.gz"
F_H5AD = BASE / "results" / "43_scDRS" / "gse190452_annotated.h5ad"
F_SCORES = BASE / "results" / "03_GSE63808_sample_LDamFerroptosis_scores.csv"

rng = np.random.default_rng(25)
notes = []

# ---------------- 1. 解析 bulk series matrix ----------------
print("[1/7] 解析 GSE63808 series matrix ...")
with gzip.open(F_SERIES, "rt", encoding="utf-8", errors="replace") as fh:
    lines_meta = []
    table_begin = None
    for i, line in enumerate(fh):
        if line.startswith("!Sample_geo_accession"):
            samples = [s.strip().strip('"') for s in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!series_matrix_table_begin"):
            table_begin = i
            break
assert table_begin is not None
expr = pd.read_csv(F_SERIES, sep="\t", skiprows=table_begin + 1, compression="gzip",
                   low_memory=False)
expr = expr[~expr.iloc[:, 0].astype(str).str.startswith("!")]
expr = expr.set_index(expr.columns[0])
expr = expr.apply(pd.to_numeric, errors="coerce")
expr.columns = [c.strip().strip('"') for c in expr.columns]
assert list(expr.columns) == samples, "样本顺序与 !Sample_geo_accession 不一致"
n_probes_raw = expr.shape[0]
log2_done = False
if np.nanmax(expr.values) > 100:
    expr = np.log2(expr + 1)
    expr = pd.DataFrame(expr, index=expr.index, columns=expr.columns)
    log2_done = True
print(f"    探针矩阵 {expr.shape[0]} x {expr.shape[1]}; log2 变换: {log2_done}; "
      f"max={np.nanmax(expr.values):.2f}")

# ---------------- 2. 探针→基因 ----------------
print("[2/7] GPL6947 探针→基因映射 ...")
with gzip.open(F_GPL, "rt", encoding="utf-8", errors="replace") as fh:
    hdr = None
    for i, line in enumerate(fh):
        if line.startswith("ID\t"):
            hdr = i
            break
gpl = pd.read_csv(F_GPL, sep="\t", skiprows=hdr, compression="gzip", low_memory=False,
                  dtype=str)
gpl = gpl[["ID", "Gene symbol"]].dropna()
gpl["symbol"] = gpl["Gene symbol"].str.split("///").str[0].str.strip()
gpl = gpl[gpl["symbol"] != ""]
probe2sym = dict(zip(gpl["ID"], gpl["symbol"]))
expr["symbol"] = expr.index.map(probe2sym)
annot = expr.dropna(subset=["symbol"])
annot = annot[annot["symbol"] != ""]
bulk_gene = annot.drop(columns=["symbol"]).groupby(annot["symbol"]).mean()
print(f"    注释成功探针 {annot.shape[0]}/{n_probes_raw}; 基因级矩阵 "
      f"{bulk_gene.shape[0]} x {bulk_gene.shape[1]}")

# ---------------- 3. 单细胞参考: 归一化 + 签名 ----------------
print("[3/7] 读取单细胞参考并构建签名矩阵 ...")
a = ad.read_h5ad(F_H5AD)
celltypes = sorted(a.obs["celltype"].astype(str).unique().tolist())
print(f"    细胞 {a.n_obs}, 细胞类型 {len(celltypes)}: {celltypes}")
if len(celltypes) < 8:
    notes.append(f"参考 h5ad 中实际只有 {len(celltypes)} 类细胞（{celltypes}）；"
                 "Pericyte 无注释细胞，签名矩阵按实际类别构建。")
an = ad.AnnData(X=a.layers["counts"].copy(), obs=a.obs[["celltype"]].copy(),
                var=a.var.copy())
sc.pp.normalize_total(an, target_sum=1e4)
sc.pp.log1p(an)

# 候选基因: sc ∩ bulk
bulk_genes = set(bulk_gene.index)
cand = [g for g in an.var_names if g in bulk_genes]
# 表达过滤: 至少在 1% 细胞中表达
Xc = an[:, cand].X
det = np.asarray((Xc > 0).sum(axis=0)).ravel() / an.n_obs
cand = [g for g, d in zip(cand, det) if d >= 0.01]
print(f"    sc∩bulk 候选基因 {len(cand)}")

# ---------------- 4. 特征基因 (top100/类, 分层子抽样 wilcoxon) ----------------
print("[4/7] rank_genes_groups 选特征基因 ...")
MAX_PER_TYPE = 3000
idx_sub = []
cts = an.obs["celltype"].astype(str).values
for ct in celltypes:
    ii = np.where(cts == ct)[0]
    if len(ii) > MAX_PER_TYPE:
        ii = rng.choice(ii, MAX_PER_TYPE, replace=False)
    idx_sub.append(ii)
idx_sub = np.sort(np.concatenate(idx_sub))
ad_sub = an[idx_sub, cand].copy()
sc.tl.rank_genes_groups(ad_sub, "celltype", method="wilcoxon", n_genes=100)
markers = []
for ct in celltypes:
    markers += list(ad_sub.uns["rank_genes_groups"]["names"][ct])
markers = sorted(set(markers))
markers = [g for g in markers if g in bulk_genes]
print(f"    合并去重特征基因 {len(markers)}（均已与 bulk 取交集）")

# 签名: 全部细胞平均 log-normalized 表达
midx = an.var_names.get_indexer(markers)
sig = np.zeros((len(markers), len(celltypes)))
for j, ct in enumerate(celltypes):
    mask = cts == ct
    col = an[mask, markers].X.mean(axis=0)
    sig[:, j] = np.asarray(col).ravel()
sig_df = pd.DataFrame(sig, index=markers, columns=celltypes)
sig_df.to_csv(RES / "signature_matrix.csv")

# ---------------- 5. NNLS 反卷积 ----------------
print("[5/7] NNLS 反卷积 (基因中心差分对数空间主分析 + 线性空间敏感性) ...")
B = bulk_gene.loc[markers].values          # log2 空间
# 基因中心差分: 签名按细胞类型均值中心化, bulk 按样本均值中心化,
# 使 NNLS 拟合基因在细胞类型间的相对差异, 消除 bulk 芯片与 sc 参考的基线/量纲差异
S_cen = sig - sig.mean(axis=1, keepdims=True)
B_cen = B - np.nanmean(B, axis=1, keepdims=True)
B_cen = np.nan_to_num(B_cen)               # 个别基因缺失按 0(=样本均值)处理
n_samples = B.shape[1]
props_log = np.zeros((n_samples, len(celltypes)))
for i in range(n_samples):
    coef, _ = nnls(S_cen, B_cen[:, i])
    ssum = coef.sum()
    props_log[i] = coef / ssum if ssum > 0 else np.nan
prop_df = pd.DataFrame(props_log, index=bulk_gene.columns, columns=celltypes)
prop_df.index.name = "sample"
prop_df.to_csv(RES / "deconv_proportions.csv")

# 敏感性: 线性空间 NNLS（expm1(均值log) 为近似; bulk 截断负值后取 2^x）
S_lin = np.expm1(np.clip(sig, 0, None))
B_lin = np.power(2.0, np.clip(B, 0, None))
props_lin = np.zeros_like(props_log)
for i in range(n_samples):
    coef, _ = nnls(S_lin, B_lin[:, i])
    ssum = coef.sum()
    props_lin[i] = coef / ssum if ssum > 0 else np.nan
prop_lin_df = pd.DataFrame(props_lin, index=bulk_gene.columns, columns=celltypes)
prop_lin_df.to_csv(RES / "deconv_proportions_linear_sensitivity.csv")

# ---------------- 6. 验证 ----------------
print("[6/7] 验证相关性 ...")
scores = pd.read_csv(F_SCORES).set_index("sample")
common = prop_df.index.intersection(scores.index)
assert len(common) == prop_df.shape[0], f"样本对齐失败: {len(common)}"
scores = scores.loc[common]
P = prop_df.loc[common]

def spr(x, y):
    r, p = spearmanr(x, y)
    return float(r), float(p)

rows = []
def add_val(name, x, y, note=""):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    n_drop = int((~ok).sum())
    if n_drop:
        note = (note + f";剔除NaN样本{n_drop}例").strip(";")
    r, p = spr(x[ok], y[ok])
    rows.append(dict(test=name, rho=r, p=p, n=int(ok.sum()), note=note))
    print(f"    {name}: rho={r:.3f}, P={p:.3g}, n={ok.sum()}")

add_val("Microglia_prop_vs_LDAM_ssGSEA", P["Microglia"], scores["LDAM"])
add_val("Microglia_prop_vs_Ferroptosis_ssGSEA", P["Microglia"], scores["Ferroptosis"])
for g in ["AIF1", "CX3CR1", "P2RY12"]:
    if g in bulk_gene.index:
        add_val(f"Microglia_prop_vs_bulk_{g}", P["Microglia"],
                bulk_gene.loc[g, common].values, note="内验证:微胶质标记基因")
    else:
        notes.append(f"bulk 矩阵中未找到 {g}")
if "SNAP25" in bulk_gene.index:
    add_val("ExcNeuron_prop_vs_bulk_SNAP25", P["ExcNeuron"],
            bulk_gene.loc["SNAP25", common].values, note="反向sanity check:突触基因")
# 敏感性: 线性空间比例 vs 对数空间比例一致性
lin_mic = prop_lin_df.loc[common, "Microglia"]
if lin_mic.nunique() <= 1:
    rows.append(dict(test="Microglia_prop_log_vs_linearNNLS", rho=np.nan, p=np.nan,
                     n=len(common),
                     note="线性空间NNLS微胶质估计退化为常数(全0), 无法相关; "
                          "该现象本身说明线性空间结果不可靠, 主分析结论以对数空间为准"))
else:
    r, p = spr(P["Microglia"], lin_mic)
    rows.append(dict(test="Microglia_prop_log_vs_linearNNLS", rho=float(r), p=float(p),
                     n=len(common), note="反卷积空间选择稳健性"))
val_df = pd.DataFrame(rows)
val_df.to_csv(RES / "deconv_validation.csv", index=False)

# ---------------- 7. 分层 ----------------
print("[7/7] 按微胶质比例 tertile 分层 ...")
mic = P["Microglia"]
t1, t2 = mic.quantile(1 / 3), mic.quantile(2 / 3)
grp = pd.Series(np.where(mic <= t1, "low", np.where(mic >= t2, "high", "mid")),
                index=mic.index)
strat_rows = []
for trait in ["LDAM", "Ferroptosis"]:
    lo = scores.loc[grp[grp == "low"].index, trait]
    hi = scores.loc[grp[grp == "high"].index, trait]
    u, p = mannwhitneyu(hi, lo, alternative="two-sided")
    strat_rows.append(dict(trait=trait, split="microglia_prop_tertile",
                           n_high=len(hi), n_low=len(lo),
                           median_high=float(hi.median()), median_low=float(lo.median()),
                           U=float(u), p=float(p)))
    print(f"    {trait}: high median={hi.median():.4f} vs low median={lo.median():.4f}, "
          f"MW P={p:.3g}")
# 敏感性: 中位数二分
med = mic.median()
grp2 = pd.Series(np.where(mic >= med, "high", "low"), index=mic.index)
for trait in ["LDAM", "Ferroptosis"]:
    lo = scores.loc[grp2[grp2 == "low"].index, trait]
    hi = scores.loc[grp2[grp2 == "high"].index, trait]
    u, p = mannwhitneyu(hi, lo, alternative="two-sided")
    strat_rows.append(dict(trait=trait, split="microglia_prop_median",
                           n_high=len(hi), n_low=len(lo),
                           median_high=float(hi.median()), median_low=float(lo.median()),
                           U=float(u), p=float(p)))
strat_df = pd.DataFrame(strat_rows)
strat_df.to_csv(RES / "stratification_results.csv", index=False)
grp.to_frame("microglia_group_tertile").to_csv(RES / "sample_stratification_groups.csv")

# ---------------- 图 ----------------
print("[图] 绘制 v25 图 ...")
import matplotlib
matplotlib.use("Agg")
try:
    sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
    from daimon_runtime import setup_plot
    setup_plot()
except Exception:
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
import matplotlib.pyplot as plt

rho_l = rows[0]["rho"]; p_l = rows[0]["p"]
rho_f = rows[1]["rho"]; p_f = rows[1]["p"]
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
for ax, trait, r_, p_, c in [(axes[0], "LDAM", rho_l, p_l, "#C0392B"),
                             (axes[1], "Ferroptosis", rho_f, p_f, "#2C5F8A")]:
    ax.scatter(P["Microglia"], scores[trait], s=12, alpha=0.6, color=c,
               edgecolors="none")
    m, b = np.polyfit(P["Microglia"], scores[trait], 1)
    xs = np.linspace(P["Microglia"].min(), P["Microglia"].max(), 50)
    ax.plot(xs, m * xs + b, color="black", lw=1, ls="--")
    ax.set_xlabel("Estimated microglia proportion (NNLS)")
    ax.set_ylabel(f"{trait} ssGSEA score (bulk)")
    ax.text(0.03, 0.97, f"Spearman ρ = {r_:.3f}\nP = {p_:.2e}",
            transform=ax.transAxes, va="top", fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "Fig_deconv_validation.png", dpi=300, bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
for k, trait in enumerate(["LDAM", "Ferroptosis"]):
    ax = axes[k]
    lo = scores.loc[grp[grp == "low"].index, trait].values
    hi = scores.loc[grp[grp == "high"].index, trait].values
    bp = ax.boxplot([lo, hi], tick_labels=["Microglia-low\n(tertile 1)",
                                           "Microglia-high\n(tertile 3)"],
                    patch_artist=True, widths=0.5, showfliers=False,
                    medianprops=dict(color="black"))
    for patch, c in zip(bp["boxes"], ["#B0C4DE", "#C0392B"]):
        patch.set_facecolor(c); patch.set_alpha(0.7)
    for xi, arr in [(1, lo), (2, hi)]:
        ax.scatter(np.full(len(arr), xi) + rng.uniform(-0.08, 0.08, len(arr)),
                   arr, s=8, alpha=0.4, color="#4d4d4d", edgecolors="none")
    p_ = strat_df[(strat_df.trait == trait) &
                  (strat_df.split == "microglia_prop_tertile")]["p"].iloc[0]
    ax.set_ylabel(f"{trait} ssGSEA score (bulk)")
    ax.set_title(f"Mann-Whitney P = {p_:.2e}", fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "Fig_deconv_stratification.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------- summary ----------------
summary = dict(
    dataset="GSE63808", platform="GPL6947 (Illumina HumanHT-12 V3)",
    n_samples=int(n_samples),
    phenotype="全部 129 例 !Sample_characteristics 均为 'phenotype: epilepsy'，"
              "无对照组，未做组间比较与生存/预后分析",
    probes_raw=int(n_probes_raw), probes_annotated=int(annot.shape[0]),
    genes_bulk=int(bulk_gene.shape[0]), log2_transformed=bool(log2_done),
    reference="results/43_scDRS/gse190452_annotated.h5ad (GSE190452 snRNA-seq)",
    celltypes=celltypes,
    n_marker_genes=len(markers),
    marker_selection="wilcoxon rank_genes_groups, 每类 top100, "
                     f"差异检验在每类≤{MAX_PER_TYPE}细胞分层子抽样上 (seed=25), "
                     "签名均值用全部细胞",
    deconv_method="基因中心差分对数空间 NNLS (scipy.optimize.nnls): 签名按细胞类型均值中心化, "
                  "bulk 按样本均值中心化, 系数归一化为比例; 线性空间 NNLS 作敏感性. "
                  "注: 绝对水平对数空间截断 NNLS 结果退化(ExcNeuron=0, Endothelial≈47%), 已弃用并记录",
    mean_proportions={ct: float(prop_df[ct].mean()) for ct in celltypes},
    validation=rows, stratification=strat_rows, notes=notes,
)
with open(RES / "deconv_summary.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)
print("完成。输出见 results/63_deconv/ 与 figures/v25/")
