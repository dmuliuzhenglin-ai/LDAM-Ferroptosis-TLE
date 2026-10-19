# -*- coding: utf-8 -*-
"""Fig_axis_cismr.png: axis-gene cis-MR screen (panels A-B).

Run:  python scripts/p28_fig_axis_cismr.py
All values are read from results/60_axis_cismr/axis_cismr_results.csv.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results" / "60_axis_cismr" / "axis_cismr_results.csv"
OUT = ROOT / "figures" / "v25" / "Fig_axis_cismr.png"

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mr = pd.read_csv(DATA)
outcomes = ["ILAE3_focalHS", "ILAE3_all", "ILAE3_focal"]
oc_label = {"ILAE3_focalHS": "Focal epilepsy + HS",
            "ILAE3_all": "All epilepsy",
            "ILAE3_focal": "Focal epilepsy (all)"}
hit_genes = ["SCD5", "PROK2"]

fig = plt.figure(figsize=(12, 4.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.30)

# Panel A: volcano for focal-HS ----------------------------------------------------
axA = fig.add_subplot(gs[0])
fh = mr[mr["outcome"] == "ILAE3_focalHS"].dropna(subset=["beta_wald", "p_wald"])
neglog = -np.log10(fh["p_wald"])
sig = fh["fdr_wald"] < 0.05
axA.scatter(fh.loc[~sig, "beta_wald"], neglog[~sig], s=6, c="#9AA5B1",
            alpha=0.5, lw=0, rasterized=True)
axA.scatter(fh.loc[sig & (fh["beta_wald"] > 0), "beta_wald"],
            neglog[sig & (fh["beta_wald"] > 0)], s=26, c="#C0392B",
            zorder=3, lw=0.4, edgecolor="black")
axA.scatter(fh.loc[sig & (fh["beta_wald"] < 0), "beta_wald"],
            neglog[sig & (fh["beta_wald"] < 0)], s=26, c="#2471A3",
            zorder=3, lw=0.4, edgecolor="black")
offsets = {"SCD5": (14, -4), "PROK2": (12, -16)}
for g in hit_genes:
    row = fh[fh["gene"] == g].iloc[0]
    axA.annotate(f"{g}\n$\\beta={row['beta_wald']:.3f}$, FDR={row['fdr_wald']:.3f}",
                 xy=(row["beta_wald"], -np.log10(row["p_wald"])),
                 xytext=offsets[g], textcoords="offset points", fontsize=8.5,
                 fontweight="bold",
                 arrowprops=dict(arrowstyle="-", lw=0.7, color="#555555"))
axA.set_ylim(top=neglog.max() * 1.32)
axA.set_xlabel("cis-MR effect size ($\\beta$, Wald ratio / IVW lead)")
axA.set_ylabel("$-\\log_{10}(P)$")
axA.set_title(f"Axis-gene cis-MR screen — {oc_label['ILAE3_focalHS']}\n"
              f"(n = {len(fh)} genes tested)", fontsize=10)
for s in ("top", "right"):
    axA.spines[s].set_visible(False)
axA.text(-0.12, 1.04, "A", transform=axA.transAxes, fontsize=14,
         fontweight="bold", va="top")

# Panel B: forest plot SCD5 / PROK2 across outcomes ----------------------------------
axB = fig.add_subplot(gs[1])
colors = {"ILAE3_focalHS": "#C0392B", "ILAE3_all": "#7F8C8D",
          "ILAE3_focal": "#2471A3"}
yticklabels, ypos, group_header_y = [], [], []
y = 0.0
for gi, g in enumerate(hit_genes):
    group_header_y.append(y - 0.55)
    for oi, oc in enumerate(outcomes):
        row = mr[(mr["gene"] == g) & (mr["outcome"] == oc)].iloc[0]
        # use IVW when available (n_inst >= 2), else lead-SNP Wald ratio
        if pd.notna(row["beta_ivw"]):
            beta, se, est = row["beta_ivw"], row["se_ivw"], "IVW"
        else:
            beta, se, est = row["beta_wald"], row["se_wald"], "Wald"
        lo, hi = beta - 1.96 * se, beta + 1.96 * se
        sig = row["fdr_wald"] < 0.05
        axB.errorbar(beta, y, xerr=[[beta - lo], [hi - beta]],
                     fmt="o", ms=6 if sig else 4.5, color=colors[oc],
                     ecolor=colors[oc], elinewidth=1.4, capsize=3,
                     mfc=colors[oc], mec="black" if sig else colors[oc],
                     mew=0.8 if sig else 0, zorder=3)
        axB.text(0.255, y,
                 f"$\\beta$={beta:+.3f} ({lo:+.3f}, {hi:+.3f}); "
                 f"FDR={row['fdr_wald']:.3f}" + (" *" if sig else ""),
                 va="center", fontsize=7.4,
                 color=colors[oc] if sig else "#555555")
        yticklabels.append(f"{oc_label[oc]}")
        ypos.append(y)
        y += 1.0
    if gi == 0:
        y += 0.9
axB.axvline(0, color="black", lw=0.8, ls="--")
axB.set_yticks(ypos)
axB.set_yticklabels(yticklabels, fontsize=8.5)
axB.set_ylim(-1.15, y - 0.4)
axB.invert_yaxis()
axB.set_xlim(-0.24, 0.72)
axB.set_xlabel("cis-MR $\\beta$ (95% CI)")
axB.set_title("SCD5 & PROK2 across ILAE3 outcomes\n(focal-HS–specific signals)",
              fontsize=10)
# gene group headers above each block
for gi, g in enumerate(hit_genes):
    axB.text(-0.225, group_header_y[gi], g, fontsize=10.5, fontweight="bold",
             va="center", color="#222222")
handles = [plt.Line2D([0], [0], marker="o", color="w", mfc=colors[o],
                      ms=6, label=oc_label[o]) for o in outcomes]
axB.legend(handles=handles, frameon=False, fontsize=7.5,
           loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
           columnspacing=0.9, handletextpad=0.4)
for s in ("top", "right"):
    axB.spines[s].set_visible(False)
axB.text(-0.14, 1.04, "B", transform=axB.transAxes, fontsize=14,
         fontweight="bold", va="top")

fig.savefig(OUT, bbox_inches="tight", dpi=300)
print("saved", OUT)
