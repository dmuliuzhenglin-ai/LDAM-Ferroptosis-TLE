# -*- coding: utf-8 -*-
"""
p30_spatial_human_figs.py
空间转录组图件绘制（Fig_spatial_human_maps / Fig_spatial_human_stats）
运行：python scripts\\p30_spatial_human_figs.py  （托管 python，300 dpi）
输入：results/65_spatial_human/ 下的 spot_scores_*.csv、stats_summary.json、
      *_tissue_lowres.png、*_scalefactors.json
"""
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "65_spatial_human"
FIG = ROOT / "figures" / "v25"
FIG.mkdir(parents=True, exist_ok=True)

SAMPLES = {"CTRL": "Control hippocampus", "HS": "MTLE-HS hippocampus"}

dfs, imgs, sfs = {}, {}, {}
for tag in SAMPLES:
    dfs[tag] = pd.read_csv(RES / f"spot_scores_{tag}.csv")
    imgs[tag] = plt.imread(RES / f"{tag}_tissue_lowres.png")
    sfs[tag] = json.load(open(RES / f"{tag}_scalefactors.json"))
stats_j = json.load(open(RES / "stats_summary.json"))

CMAP = "viridis"


def spot_xy(tag):
    """fullres 像素坐标 -> lowres 图像坐标。"""
    df = dfs[tag]
    sf = sfs[tag]["tissue_lowres_scalef"]
    return df["pxl_col_in_fullres"] * sf, df["pxl_row_in_fullres"] * sf


def spatial_ax(ax, tag, values, title, vmin=None, vmax=None):
    """在 lowres 组织图上叠加 spot 散点。"""
    x, y = spot_xy(tag)
    ax.imshow(imgs[tag])
    sc = ax.scatter(x, y, c=values, cmap=CMAP, s=6, vmin=vmin, vmax=vmax,
                    linewidths=0, alpha=0.95)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    return sc


# ============ Fig_spatial_human_maps ============
# A/B: LDAM; C/D: Ferroptosis; E/F: LGALS3 表达（两样本共用同一色标范围）
fig, axes = plt.subplots(3, 2, figsize=(9, 12))
panels = [("score_LDAM", "LDAM score"),
          ("score_Ferroptosis", "Ferroptosis score"),
          ("expr_LGALS3", "LGALS3 expression (log1p)")]
letters = "ABCDEF"
for i, (col, label) in enumerate(panels):
    lo = min(dfs[t][col].quantile(0.02) for t in SAMPLES)
    hi = max(dfs[t][col].quantile(0.98) for t in SAMPLES)
    for j, tag in enumerate(SAMPLES):
        ax = axes[i, j]
        im = spatial_ax(ax, tag, dfs[tag][col], f"{SAMPLES[tag]}\n{label}",
                        vmin=lo, vmax=hi)
        ax.text(-0.05, 1.03, letters[2 * i + j], transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="top")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=7)
fig.suptitle("Spatial maps: LDAM / ferroptosis axis in human hippocampus "
             "(Visium, proof-of-concept n=1 vs n=1)", fontsize=12, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.99])
fig.savefig(FIG / "Fig_spatial_human_maps.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("[fig] Fig_spatial_human_maps.png")

# ============ Fig_spatial_human_stats ============
fig = plt.figure(figsize=(15, 4.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1.0, 1.2], wspace=0.3)

# --- Panel A: 轴分数分布小提琴图 ---
axA = fig.add_subplot(gs[0])
axes_names = ["score_LDAM", "score_Ferroptosis", "score_Module5"]
labels = ["LDAM", "Ferroptosis", "Module5"]
colors = {"CTRL": "#4C9BD6", "HS": "#D6604D"}
pos = np.arange(len(axes_names))
for k, tag in enumerate(SAMPLES):
    data = [dfs[tag][c] for c in axes_names]
    vp = axA.violinplot(data, positions=pos + (k - 0.5) * 0.36, widths=0.32,
                        showmedians=True, showextrema=False)
    for b in vp["bodies"]:
        b.set_facecolor(colors[tag]); b.set_alpha(0.7)
    vp["cmedians"].set_color("black")
    axA.scatter([], [], color=colors[tag], s=40, label=tag)
for i, c in enumerate(axes_names):
    p = stats_j["comparisons"][c]["p_value"]
    r = stats_j["comparisons"][c]["rank_biserial_r"]
    ymax = max(dfs[t][c].max() for t in SAMPLES)
    ptxt = "P<1e-300" if p == 0 else f"P={p:.1e}"
    axA.text(i, ymax * 1.04, f"{ptxt}\nr={r:+.2f}", ha="center",
             va="bottom", fontsize=7)
axA.set_xticks(pos); axA.set_xticklabels(labels)
axA.set_ylabel("Gene-set score")
axA.set_ylim(top=max(dfs[t][c].max() for t in SAMPLES
                     for c in axes_names) * 1.3)
axA.legend(frameon=False, fontsize=8, loc="upper left")
axA.set_title("A  Axis score distributions (Mann-Whitney)",
              fontsize=10, loc="left", fontweight="bold")

# --- Panel B: LDAM 分数 vs 微胶质富集度散点 ---
axB = fig.add_subplot(gs[1])
for tag in SAMPLES:
    df = dfs[tag]
    axB.scatter(df["microglia_enrichment"], df["score_LDAM"], s=3,
                alpha=0.35, color=colors[tag], label=tag, rasterized=True)
    rho = stats_j["spearman"][f"score_LDAM__{tag}"]["rho"]
    p = stats_j["spearman"][f"score_LDAM__{tag}"]["p_value"]
    # 线性拟合辅助线
    m, b = np.polyfit(df["microglia_enrichment"], df["score_LDAM"], 1)
    xs = np.linspace(df["microglia_enrichment"].min(),
                     df["microglia_enrichment"].max(), 50)
    axB.plot(xs, m * xs + b, color=colors[tag], lw=1.5)
    print(f"[spearman] {tag}: rho={rho:.3f}, P={p:.1e}")
axB.text(0.03, 0.97,
         "\n".join(f"{t}: ρ={stats_j['spearman'][f'score_LDAM__{t}']['rho']:+.3f}, "
                   f"P={stats_j['spearman'][f'score_LDAM__{t}']['p_value']:.1e}"
                   for t in SAMPLES),
         transform=axB.transAxes, fontsize=8, va="top")
axB.set_xlabel("Microglia enrichment (mean AIF1/CX3CR1/P2RY12/TMEM119)")
axB.set_ylabel("LDAM score")
axB.legend(frameon=False, fontsize=8, loc="lower right")
axB.set_title("B  LDAM score vs microglia enrichment (Spearman)",
              fontsize=10, loc="left", fontweight="bold")

# --- Panel C: 各海马亚区轴分数条形图 ---
axC = fig.add_subplot(gs[2])
regions = ["DG-like (PROX1+)", "CA1-like (NETO1+)", "astro-enriched",
           "microglia-enriched", "mixed/other"]
disp = ["DG-like", "CA1-like", "Astro-enriched", "Microglia-enriched",
        "Mixed/other"]
x = np.arange(len(regions))
w = 0.38
for k, tag in enumerate(SAMPLES):
    rs = stats_j["region_summary"][tag]
    vals = [rs[r]["LDAM_mean"] if r in rs else np.nan for r in regions]
    axC.bar(x + (k - 0.5) * w, vals, w, color=colors[tag], alpha=0.85,
            label=tag)
axC.set_xticks(x)
axC.set_xticklabels(disp, rotation=20, ha="right", fontsize=8)
axC.set_ylabel("Mean LDAM score")
axC.legend(frameon=False, fontsize=8)
axC.set_title("C  LDAM score by hippocampal subregion\n"
              "(expression-cluster annotation, heuristic)",
              fontsize=10, loc="left", fontweight="bold")

fig.tight_layout()
fig.savefig(FIG / "Fig_spatial_human_stats.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("[fig] Fig_spatial_human_stats.png")
print("DONE")
