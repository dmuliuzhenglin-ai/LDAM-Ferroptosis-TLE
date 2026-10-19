"""v2.4 figure: E16.5 embryo gsMap Cauchy results for ILAE2018_all and ILAE3 focal-HS."""
import sys, os
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from daimon_runtime import setup_plot
setup_plot()

BASE = os.getcwd()
FIG = os.path.join(BASE, "figures", "replication_v24")

def load(tag):
    sub = "quick_mode_E16_ILAE2018_all" if tag == "ILAE2018_all" else "quick_mode_E16_focal_HS"
    f = os.path.join(BASE, "results", "53_ILAE2023_gsMap", sub, "E16.5_E1S1",
                     "cauchy_combination", f"E16.5_E1S1_{tag}.Cauchy.csv.gz")
    return pd.read_csv(f, compression="gzip")

a = load("ILAE2018_all").rename(columns={"p_cauchy": "p_2018"})[["annotation", "p_2018"]]
b = load("focal_HS").rename(columns={"p_cauchy": "p_focalHS"})[["annotation", "p_focalHS"]]
d = a.merge(b, on="annotation")
d = d.sort_values("p_2018", ascending=False).reset_index(drop=True)

y = np.arange(len(d))
fig, ax = plt.subplots(figsize=(7.2, 6.2))
ax.scatter(-np.log10(d.p_2018), y, s=22, c="#D64541", label="ILAE 2018 all epilepsy", zorder=3)
ax.scatter(-np.log10(d.p_focalHS), y, s=22, c="#4C9BD6", label="ILAE 2023 focal-HS", zorder=3)
thr = -np.log10(0.05 / len(d))
ax.axvline(-np.log10(0.05), c="k", lw=0.7, ls="--")
ax.text(-np.log10(0.05) + 0.03, len(d) - 0.5, "P = 0.05", fontsize=7, rotation=90, va="top")
ax.axvline(thr, c="k", lw=0.7, ls=":")
ax.text(thr + 0.03, len(d) - 0.5, "Bonferroni", fontsize=7, rotation=90, va="top")
ax.set_yticks(y, d.annotation, fontsize=8)
ax.set_xlabel(r"$-\log_{10}$ Cauchy P (organ-level enrichment)")
ax.set_title("E16.5 mouse embryo (gsMap, 121,767 spots)", fontsize=10)
ax.legend(fontsize=8, frameon=False, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "Figure_E16_gsmap_v24.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIG, "Figure_E16_gsmap_v24.pdf"), bbox_inches="tight")
print("saved")
