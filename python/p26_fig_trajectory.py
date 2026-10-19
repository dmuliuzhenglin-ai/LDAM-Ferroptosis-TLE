# -*- coding: utf-8 -*-
"""Fig_trajectory.png: microglia DPT pseudotime trajectory (panels A-C).

Run with the managed python:  python scripts/p26_fig_trajectory.py
Data export from the h5ad is delegated to gsmap_env python via subprocess.
All annotation numbers are read from trajectory_summary.json / the h5ad itself.
"""
import sys, json, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results" / "58_trajectory"
OUT = ROOT / "figures" / "v25" / "Fig_trajectory.png"
EXPORT = ROOT / "figures" / "v25" / "_traj_export.csv"
GSMAP_PY = ROOT / "gsmap_env" / "Scripts" / "python.exe"

# ---------------------------------------------------------------- export ----
EXPORT_CODE = r"""
import anndata as ad
import pandas as pd
from pathlib import Path
root = Path(sys.argv[1])
a = ad.read_h5ad(root / 'results/58_trajectory/microglia_dpt.h5ad')
df = pd.DataFrame({
    'umap1': a.obsm['X_umap'][:, 0],
    'umap2': a.obsm['X_umap'][:, 1],
    'dpt': a.obs['dpt_pseudotime'].to_numpy(),
    'LDAM': a.obs['LDAM_score'].to_numpy(),
    'ferro': a.obs['Ferroptosis_score'].to_numpy(),
    'condition': a.obs['condition'].astype(str).to_numpy(),
})
df.to_csv(root / 'figures/v25/_traj_export.csv', index=False)
print('exported', df.shape)
"""
if not EXPORT.exists():
    subprocess.run([str(GSMAP_PY), "-c",
                    "import sys; " + EXPORT_CODE.replace("sys.argv[1]", "sys.argv[1]"),
                    str(ROOT)], check=True)

# ------------------------------------------------------------------ plot ----
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

df = pd.read_csv(EXPORT)
summ = json.loads((DATA / "trajectory_summary.json").read_text(encoding="utf-8"))
rho_l, p_l = summ["spearman_dpt_LDAM"]
rho_f, p_f = summ["spearman_dpt_ferroptosis"]
n_cells = summ["n_cells"]
late_htle = summ["condition_x_tercile"]["HTLE"]["late"]


def lowess(x, y, frac=0.3):
    """Locally weighted linear regression (tricube kernel), numpy only."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    order = np.argsort(x); x, y = x[order], y[order]
    n = len(x); k = max(int(frac * n), 3)
    xs = np.linspace(x.min(), x.max(), 200)
    ys = np.empty_like(xs)
    for i, xi in enumerate(xs):
        d = np.abs(x - xi)
        h = np.partition(d, k - 1)[k - 1]
        w = np.where(d < h, (1 - (d / h) ** 3) ** 3, 0.0)
        X = np.column_stack([np.ones(n), x - xi])
        W = w[:, None]
        beta = np.linalg.lstsq(X * W, y * w, rcond=None)[0]
        ys[i] = beta[0]
    return xs, ys


def fmt_p(p):
    m, e = f"{p:.1e}".split("e")
    return rf"$P={m}\times10^{{{int(e)}}}$"


fig = plt.figure(figsize=(14, 4.3))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.35], wspace=0.52)

# Panel A: UMAP coloured by DPT ------------------------------------------------
axA = fig.add_subplot(gs[0])
rng = np.random.default_rng(0)
idx = rng.permutation(len(df))
sc = axA.scatter(df["umap1"].to_numpy()[idx], df["umap2"].to_numpy()[idx],
                 c=df["dpt"].to_numpy()[idx], cmap="viridis", s=3, lw=0,
                 rasterized=True)
axA.set_xlabel("UMAP 1"); axA.set_ylabel("UMAP 2")
axA.set_title("DPT pseudotime", fontsize=10)
cb = fig.colorbar(sc, ax=axA, fraction=0.040, pad=0.06)
cb.set_label("DPT pseudotime", fontsize=8); cb.ax.tick_params(labelsize=7)
axA.text(-0.10, 1.06, "A", transform=axA.transAxes, fontsize=14,
         fontweight="bold", va="top")

# Panel B: UMAP coloured by condition -------------------------------------------
axB = fig.add_subplot(gs[1])
colors = {"Control": "#4C72B0", "HTLE": "#C44E52"}
for cond in ["Control", "HTLE"]:
    m = (df["condition"] == cond).to_numpy()
    axB.scatter(df["umap1"].to_numpy()[m], df["umap2"].to_numpy()[m],
                c=colors[cond], s=3, lw=0, alpha=0.75, label=cond,
                rasterized=True)
axB.set_xlabel("UMAP 1"); axB.set_ylabel("UMAP 2")
axB.set_title(f"Condition (late tercile: {late_htle*100:.1f}% HTLE)",
              fontsize=10)
axB.legend(frameon=False, markerscale=4, fontsize=8, loc="best")
axB.text(-0.10, 1.06, "B", transform=axB.transAxes, fontsize=14,
         fontweight="bold", va="top")

# Panel C: scores along DPT ------------------------------------------------------
axC = fig.add_subplot(gs[2])
x = df["dpt"].to_numpy()
axC.scatter(x, df["LDAM"], s=3, c="#C9A227", alpha=0.25, lw=0, rasterized=True)
xs, ys = lowess(x, df["LDAM"].to_numpy())
axC.plot(xs, ys, color="#8C6D1F", lw=2.2,
         label=rf"LDAM score  $\rho={rho_l:.3f}$, " + fmt_p(p_l))
axC.scatter(x, df["ferro"], s=3, c="#55A3A8", alpha=0.25, lw=0, rasterized=True)
xs, ys = lowess(x, df["ferro"].to_numpy())
axC.plot(xs, ys, color="#1F6F74", lw=2.2,
         label=rf"Ferroptosis score  $\rho={rho_f:.3f}$, " + fmt_p(p_f))
axC.set_xlabel("DPT pseudotime")
axC.set_ylabel("Module score")
axC.set_title(f"Axis scores along pseudotime (n = {n_cells} microglia)",
              fontsize=10)
axC.legend(frameon=False, fontsize=8, loc="upper left")
axC.text(-0.08, 1.06, "C", transform=axC.transAxes, fontsize=14,
         fontweight="bold", va="top")

for ax in (axA, axB, axC):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

fig.savefig(OUT, bbox_inches="tight", dpi=300)
print("saved", OUT)
