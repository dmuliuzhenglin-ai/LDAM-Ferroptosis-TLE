# -*- coding: utf-8 -*-
"""p23_docking_plot.py — bar chart of Vina affinities -> figures/v25/Fig_docking_affinity.png
Run with managed python (has matplotlib + daimon_runtime)."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "results" / "62_docking" / "docking_affinities.csv"
FIG = ROOT / "figures" / "v25" / "Fig_docking_affinity.png"
FIG.parent.mkdir(parents=True, exist_ok=True)

GROUP = {  # ligand -> (display label, color)
    "TD139":        ("TD139 (positive control)", "#C0392B"),
    "lactose":      ("Lactose (natural ligand)", "#7F8C8D"),
    "daunorubicin": ("Daunorubicin", "#2E86C1"),
    "dorsomorphin": ("Dorsomorphin", "#2E86C1"),
    "mitoxantrone": ("Mitoxantrone", "#2E86C1"),
    "triptolide":   ("Triptolide", "#2E86C1"),
}

rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
data = [(r["ligand"], float(r["best_affinity_kcal_mol"])) for r in rows]
data.sort(key=lambda x: x[1])  # most negative first -> top of chart

labels = [GROUP[k][0] for k, _ in data]
vals = [v for _, v in data]
colors = [GROUP[k][1] for k, _ in data]

fig, ax = plt.subplots(figsize=(7.2, 4.2))
bars = ax.barh(range(len(vals)), vals, color=colors, edgecolor="black", linewidth=0.6, height=0.62)
ax.set_yticks(range(len(vals)))
ax.set_yticklabels(labels, fontsize=10)
ax.invert_yaxis()  # best (lowest) on top
for i, v in enumerate(vals):
    ax.text(v - 0.06, i, f"{v:.2f}", va="center", ha="right", fontsize=9.5, fontweight="bold")
xmin = min(vals)
ax.set_xlim(xmin - 0.9, 0)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("Vina binding affinity (kcal/mol, lower = stronger)", fontsize=11)
ax.set_title("Molecular docking of candidate compounds to Galectin-3 (PDB 1KJL)", fontsize=11.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="x", labelsize=10)
fig.tight_layout()
fig.savefig(FIG, dpi=300, bbox_inches="tight")
print(f"[fig] -> {FIG}")
