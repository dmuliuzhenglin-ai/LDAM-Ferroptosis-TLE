"""Forest plot: microglia (Bryois 2022) cell-type-specific cis-MR, focal-HS outcome.
All 12 axis genes; untested-in-microglia genes shown as grey annotation rows.
Colour: nominal P<0.05 vs not significant (all null here -> shown honestly).
Output: figures/v25/Fig_microglia_MR.png (300 dpi)
"""
import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

BASE = os.getcwd()
OUT = os.path.join(BASE, "figures", "v25")
os.makedirs(OUT, exist_ok=True)

R = pd.read_csv(os.path.join(BASE, "results", "64_microglia_cismr", "microglia_cismr_results.csv"))
fhs = R[R.outcome == "ILAE3_focalHS"].set_index("gene")
comp = pd.read_csv(os.path.join(BASE, "results", "64_microglia_cismr", "comparison_with_blood.csv"))
blood = comp[comp.outcome == "ILAE3_focalHS"].set_index("gene")

GENE_ORDER = ["S100A10", "SYNGR2", "PROK2", "SCD5", "IKZF1", "JAZF1",
              "CD9", "ETV6", "LGALS3", "IRF8", "LPCAT2", "RUNX1"]

fig, ax = plt.subplots(figsize=(9.5, 6.5))
ypos = np.arange(len(GENE_ORDER))[::-1]
for y, g in zip(ypos, GENE_ORDER):
    if g in fhs.index and pd.notna(fhs.loc[g, "beta_mr"]):
        r = fhs.loc[g]
        b, se, p = r.beta_mr, r.se_mr, r.p_mr
        lo, hi = b - 1.96 * se, b + 1.96 * se
        sig = p < 0.05
        col = "#C0392B" if sig else "#2C3E50"
        ax.plot([lo, hi], [y, y], color=col, lw=2.2, solid_capstyle="round", zorder=3)
        ax.scatter([b], [y], s=90 + 4 * r.F, color=col, zorder=4,
                   marker="D" if sig else "o", edgecolor="white", linewidth=0.8)
        lab = f"P={p:.3f}"
        tier = "GW-sig IV" if "genome-wide" in r.instrument_tier else \
               "suggestive IV" if "suggestive" in r.instrument_tier else "sub-threshold IV"
        ax.text(1.02, y, f"{r.lead_rsid}  F={r.F:.1f} ({tier})  {lab}",
                transform=ax.get_yaxis_transform(), fontsize=8, va="center", color=col)
    else:
        note = "not tested in microglia panel"
        if g in blood.index and pd.notna(blood.loc[g, "blood_p_wald"]):
            br = blood.loc[g]
            note += f"  |  blood eQTLGen: beta={br.blood_beta_wald:.3f}, P={br.blood_p_wald:.3g}"
            if br.blood_fdr_wald < 0.05:
                note += " (FDR<0.05 in blood)"
        ax.text(0.02, y, note, transform=ax.get_yaxis_transform(),
                fontsize=8, va="center", style="italic", color="#7F8C8D",
                bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.2))
        ax.scatter([0], [y], s=40, marker="x", color="#BDC3C7", zorder=3)

ax.axvline(0, color="#95A5A6", lw=1, ls="--")
ax.set_yticks(ypos)
ax.set_yticklabels(GENE_ORDER, fontsize=10)
ax.set_xlabel("Microglia cis-MR beta on focal epilepsy with hippocampal sclerosis (Wald ratio, 95% CI)", fontsize=10)
ax.set_title("Cell-type-specific cis-MR: Bryois 2022 microglia eQTL (n=192) -> ILAE 2023 focal-HS\n"
             "lead-SNP Wald ratio; 6/12 axis genes absent from the microglia eQTL panel",
             fontsize=11)
xmax = max(0.5, np.nanmax(np.abs(fhs.beta_mr + 2 * fhs.se_mr)) * 1.15)
ax.set_xlim(-xmax, xmax)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "Fig_microglia_MR.png"), dpi=300, bbox_inches="tight")
print("saved", os.path.join(OUT, "Fig_microglia_MR.png"))
