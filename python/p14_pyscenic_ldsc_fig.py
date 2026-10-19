"""v2.4 figures: pySCENIC regulon scatter + LDSC rg forest plot."""
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
os.makedirs(FIG, exist_ok=True)

# ---------- Figure: pySCENIC differential regulons ----------
d = pd.read_csv(os.path.join(BASE, "results", "56_pySCENIC", "regulon_differential.csv"))
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

ax = axes[0]
sig_both = (d.fdr_ldam < 0.05) & (d.fdr_htle < 0.05)
sig_any = (d.fdr_ldam < 0.05) | (d.fdr_htle < 0.05)
ax.scatter(d.loc[~sig_any, "delta_htle"], d.loc[~sig_any, "delta_ldam"],
           s=8, c="#bbbbbb", alpha=0.6, linewidths=0, label="NS")
ax.scatter(d.loc[sig_any & ~sig_both, "delta_htle"], d.loc[sig_any & ~sig_both, "delta_ldam"],
           s=12, c="#4C9BD6", alpha=0.8, linewidths=0, label="FDR<0.05 in one contrast")
ax.scatter(d.loc[sig_both, "delta_htle"], d.loc[sig_both, "delta_ldam"],
           s=16, c="#D64541", alpha=0.9, linewidths=0, label="FDR<0.05 in both")
for _, r in d.head(0).iterrows():
    pass
top = d.sort_values("fdr_ldam").head(6)
off = {"IKZF1+": (6, 4), "ETV6+": (6, -2), "JAZF1+": (6, 2), "RUNX1+": (6, -6),
       "FLI1+": (6, 3), "FOXP2+": (6, -3)}
for _, r in top.iterrows():
    dx, dy = off.get(r["regulon"], (5, 3))
    ax.annotate(r["regulon"].replace("+", "(+)"), (r["delta_htle"], r["delta_ldam"]),
                textcoords="offset points", xytext=(dx, dy), fontsize=7.5)
irf8 = d[d.regulon == "IRF8+"].iloc[0]
ax.annotate("IRF8(+)", (irf8["delta_htle"], irf8["delta_ldam"]),
            textcoords="offset points", xytext=(6, 4), fontsize=7.5)
ax.axhline(0, c="k", lw=0.5); ax.axvline(0, c="k", lw=0.5)
ax.set_xlabel("ΔAUCell (HTLE vs control microglia)")
ax.set_ylabel("ΔAUCell (LDAM-high vs LDAM-low)")
ax.set_title("a  Regulon activity shifts", loc="left", fontweight="bold")
ax.legend(fontsize=7, frameon=False, loc="lower right")

# ---------- LDSC rg forest ----------
ax = axes[1]
rg = pd.read_csv(os.path.join(BASE, "results", "57_LDSC_rg", "ldsc_rg_results.csv"))
ext = rg[rg.trait2.isin(["AD", "STROKE", "LDL", "CRP"])].copy()
order = [("ILAE3_all", "STROKE"), ("ILAE3_all", "CRP"), ("ILAE3_all", "LDL"), ("ILAE3_all", "AD"),
         ("ILAE2018_all", "STROKE"), ("ILAE2018_all", "CRP"), ("ILAE2018_all", "LDL"), ("ILAE2018_all", "AD"),
         ("ILAE3_focal", "STROKE"), ("ILAE3_focal", "CRP"), ("ILAE3_focal", "LDL"), ("ILAE3_focal", "AD"),
         ("ILAE3_focalHS", "STROKE"), ("ILAE3_focalHS", "CRP"), ("ILAE3_focalHS", "LDL"), ("ILAE3_focalHS", "AD")]
labmap = {"ILAE3_all": "ILAE 2023 all", "ILAE2018_all": "ILAE 2018 all",
          "ILAE3_focal": "ILAE 2023 focal", "ILAE3_focalHS": "ILAE 2023 focal-HS",
          "STROKE": "Stroke", "CRP": "CRP", "LDL": "LDL-C", "AD": "Alzheimer's"}
rows = []
for t1, t2 in order:
    r = ext[(ext.trait1 == t1) & (ext.trait2 == t2)]
    if len(r):
        rows.append((f"{labmap[t1]}  vs  {labmap[t2]}", r.iloc[0].rg, r.iloc[0].se, r.iloc[0].fdr))
rows = rows[::-1]
y = np.arange(len(rows))
for i, (lab, g, se, fdr) in enumerate(rows):
    c = "#D64541" if fdr < 0.05 else "#555555"
    ax.errorbar(g, i, xerr=1.96 * se, fmt="o", ms=4, color=c, ecolor=c, elinewidth=1, capsize=2)
ax.set_yticks(y, [r[0] for r in rows], fontsize=7.5)
ax.axvline(0, c="k", lw=0.6, ls="--")
ax.set_xlabel("Genetic correlation $r_g$ (95% CI)")
ax.set_title("b  Cross-trait LDSC", loc="left", fontweight="bold")

fig.tight_layout()
fig.savefig(os.path.join(FIG, "Figure_pyscenic_ldsc_v24.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIG, "Figure_pyscenic_ldsc_v24.pdf"), bbox_inches="tight")
print("saved", os.path.join(FIG, "Figure_pyscenic_ldsc_v24.png"))
