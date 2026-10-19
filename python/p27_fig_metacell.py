# -*- coding: utf-8 -*-
"""Fig_metacell_modules.png: metacell co-expression modules (panels A-C).

Run:  python scripts/p27_fig_metacell.py
All values are read from results/59_metacell_wgcna/*.csv.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results" / "59_metacell_wgcna"
OUT = ROOT / "figures" / "v25" / "Fig_metacell_modules.png"

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

mta = pd.read_csv(DATA / "module_trait_associations.csv")
overlay = pd.read_csv(DATA / "module_regulon_overlay.csv")
members = pd.read_csv(DATA / "metacell_modules.csv")

modules = sorted(mta["module"].unique())
traits = ["LDAM", "Ferroptosis", "HTLE_vs_Control"]
trait_label = {"LDAM": "LDAM score\n(Spearman $\\rho$)",
               "Ferroptosis": "Ferroptosis score\n(Spearman $\\rho$)",
               "HTLE_vs_Control": "HTLE vs Control\n(t-statistic)"}

def stars(fdr):
    return "***" if fdr < 1e-3 else "**" if fdr < 1e-2 else "*" if fdr < 0.05 else ""

# wide tables
val = pd.DataFrame(index=modules, columns=traits, dtype=float)
ann = pd.DataFrame(index=modules, columns=traits, dtype=object)
for _, r in mta.iterrows():
    val.loc[r["module"], r["trait"]] = r["rho"]
    ann.loc[r["module"], r["trait"]] = f"{r['rho']:.2f}{stars(r['fdr'])}"

fig = plt.figure(figsize=(13.5, 4.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.05, 1.0], wspace=0.42)

# Panel A: module x trait heatmap (per-column diverging scale) -------------------
axA = fig.add_subplot(gs[0])
M = val.to_numpy(float).astype(float)
colnorm = [TwoSlopeNorm(vmin=-np.nanmax(np.abs(M[:, j])), vcenter=0,
                        vmax=np.nanmax(np.abs(M[:, j])))
           for j in range(M.shape[1])]
Rgba = plt.get_cmap("RdBu_r")
img = np.zeros((*M.shape, 4))
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        img[i, j] = Rgba(colnorm[j](M[i, j]))
axA.imshow(img, aspect="auto")
axA.grid(False)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        rgba = img[i, j]
        lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
        axA.text(j, i, ann.iloc[i, j], ha="center", va="center", fontsize=8,
                 color="white" if lum < 0.55 else "black")
axA.set_xticks(range(len(traits)))
axA.set_xticklabels([trait_label[t] for t in traits], fontsize=8)
axA.set_yticks(range(len(modules)))
axA.set_yticklabels([f"M{m} (n={int(mta[mta['module']==m]['n_genes'].iloc[0])})"
                     for m in modules], fontsize=8)
axA.set_title("Module–trait associations", fontsize=10)
axA.tick_params(length=0)
for s in axA.spines.values():
    s.set_visible(False)
axA.set_xticks(np.arange(-.5, len(traits), 1), minor=True)
axA.set_yticks(np.arange(-.5, len(modules), 1), minor=True)
axA.grid(which="minor", color="white", lw=1.5)
axA.text(-0.28, 1.06, "A", transform=axA.transAxes, fontsize=14,
         fontweight="bold", va="top")

# Panel B: module 5 member genes ---------------------------------------------------
axB = fig.add_subplot(gs[1])
m5 = sorted(members[members["module"] == 5]["gene"].tolist())
m5_rho = mta[(mta["module"] == 5) & (mta["trait"] == "LDAM")].iloc[0]
key_genes = ["APOE", "C1QB", "CD74", "HLA-DRA", "FTH1", "FTL", "CST3", "B2M"]
axB.axis("off")
axB.set_title(f"Module 5 member genes (n = {len(m5)})\n"
              rf"LDAM $\rho={m5_rho['rho']:.2f}$, "
              f"FDR={m5_rho['fdr']:.1e}", fontsize=10, loc="left")
# render gene names in columns, key genes highlighted
ncol = 4
for k, g in enumerate(m5):
    r, c = divmod(k, ncol)
    x = 0.02 + c * 0.26
    y = 0.80 - r * 0.105
    if g in key_genes:
        axB.text(x, y, g, transform=axB.transAxes, fontsize=9.5,
                 fontweight="bold", color="#B03A2E")
    else:
        axB.text(x, y, g, transform=axB.transAxes, fontsize=9,
                 color="#333333")
axB.text(0.02, 0.80 - ((len(m5) + ncol - 1) // ncol) * 0.105 - 0.04,
         "Bold red: lipid-handling / complement /\niron-storage / antigen-presentation genes",
         transform=axB.transAxes, fontsize=7.5, color="#666666")
axB.text(-0.10, 1.06, "B", transform=axB.transAxes, fontsize=14,
         fontweight="bold", va="top")

# Panel C: IRF8(+) regulon overlay per module ---------------------------------------
axC = fig.add_subplot(gs[2])
irf8 = overlay[overlay["regulon"] == "IRF8(+)"].set_index("module").loc[modules]
cols = ["#B03A2E" if m == 5 else "#8FA6BF" for m in modules]
bars = axC.bar(range(len(modules)), irf8["rho"], color=cols, width=0.65)
for i, (m, r) in enumerate(irf8.iterrows()):
    if m == 5:
        continue  # annotated separately
    y = r["rho"] + (0.03 if r["rho"] >= 0 else -0.09)
    axC.text(i, y, f"{r['rho']:.2f}", ha="center", fontsize=8)
axC.axhline(0, color="black", lw=0.8)
axC.set_xticks(range(len(modules)))
axC.set_xticklabels([f"M{m}" for m in modules], fontsize=9)
axC.set_ylim(-0.80, 0.95)
axC.set_ylabel("Spearman $\\rho$ with IRF8(+) regulon")
axC.set_title("IRF8(+) regulon activity overlay", fontsize=10)
m5_row = irf8.loc[5]
axC.annotate(rf"M5: $\rho={m5_row['rho']:.3f}$, $P={m5_row['p']:.1e}$",
             xy=(modules.index(5), m5_row["rho"] + 0.02),
             xytext=(modules.index(5) + 0.75, 0.86),
             fontsize=8.5, color="#B03A2E", ha="left",
             arrowprops=dict(arrowstyle="-", color="#B03A2E", lw=0.8))
for s in ("top", "right"):
    axC.spines[s].set_visible(False)
axC.text(-0.16, 1.06, "C", transform=axC.transAxes, fontsize=14,
         fontweight="bold", va="top")

fig.savefig(OUT, bbox_inches="tight", dpi=300)
print("saved", OUT)
