"""Regenerate panel_B.png including cortex eQTL MR rows for all 5 genes."""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
RES = os.path.join(BASE, "results", "42_MR_coloc")
FIG = os.path.join(BASE, "figures", "main_figures_v2", "fig9_panels")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "figure.facecolor": "white", "axes.facecolor": "white"})

rows = []
def add(label, mrkey_dict, relaxed=False):
    mr = mrkey_dict.get("MR:IVW (multiplicative RE)") or mrkey_dict.get("MR:Wald ratio")
    if mr:
        rows.append(dict(label=label + (" *" if relaxed else ""),
                         b=mr["b"], se=mr["se"], p=mr["p"]))

lg_all = {r["tissue"]: r for r in json.load(open(os.path.join(RES, "LGALS3_results.json")))}
lg_foc = {r["tissue"]: r for r in json.load(open(os.path.join(RES, "LGALS3_focal_results.json")))}
hub_ctx = json.load(open(os.path.join(RES, "hub_genes_cortex_results.json")))

r = lg_all["MetaBrain cortex EUR (N=2693)"]
add(f"LGALS3 eQTL cortex -> all epilepsy (nIV={r['n_instruments']})", r)
r = lg_foc["MetaBrain cortex EUR (N=2693)"]
add(f"LGALS3 eQTL cortex -> focal epilepsy (nIV={r['n_instruments']})", r)
for g in hub_ctx:
    add(f"{g['gene']} eQTL cortex -> all epilepsy (nIV={g['n_instruments']})", g,
        relaxed="RELAXED P<1e-5" in g["instrument_mode"])

for f, tag in [("15_MR_LGALS3_EPI_ALL_vs_ieu-b-8.csv", "LGALS3 pQTL (deCODE) -> all epilepsy"),
               ("15_MR_LGALS3_EPI_FOCAL_vs_ieu-b-10.csv", "LGALS3 pQTL (deCODE) -> focal epilepsy"),
               ("15_MR_LGALS3_EPI_HS_vs_ieu-b-14.csv", "LGALS3 pQTL (deCODE) -> focal + HS")]:
    d = pd.read_csv(os.path.join(BASE, "results", f))
    d = d[d["method"] == "Inverse variance weighted"].iloc[0]
    rows.append(dict(label=f"{tag} (nIV={int(d['nsnp'])})", b=d["b"], se=d["se"], p=d["pval"]))

fr = pd.DataFrame(rows)
fr["OR"] = np.exp(fr["b"]); fr["lo"] = np.exp(fr["b"] - 1.96 * fr["se"]); fr["hi"] = np.exp(fr["b"] + 1.96 * fr["se"])
fig, ax = plt.subplots(figsize=(8.8, 4.6))
yy = np.arange(len(fr))[::-1]
ncol = {True: "#8856a7", False: "#2166ac"}
for yi, (i, r) in zip(yy, fr.iterrows()):
    c = "#8856a7" if "pQTL" in r["label"] else "#2166ac"
    ax.plot([r["lo"], r["hi"]], [yi, yi], color=c, lw=1.6, solid_capstyle="round")
    ax.plot(r["OR"], yi, "s", color=c, ms=7)
    ax.text(1.16, yi, f"OR={r['OR']:.3f} ({r['lo']:.3f}-{r['hi']:.3f}), P={r['p']:.3g}",
            va="center", fontsize=8.5)
ax.axvline(1, ls="--", color="#b2182b", lw=0.9)
ax.set_yticks(yy); ax.set_yticklabels(fr["label"], fontsize=8.5)
ax.set_xscale("log"); ax.set_xlim(0.88, 1.15)
ax.set_xticks([0.9, 0.95, 1.0, 1.05, 1.1, 1.15])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel("Odds ratio for epilepsy per unit higher exposure (log scale)")
ax.set_title("Mendelian randomization: hub genes / LGALS3 (eQTL, pQTL) -> epilepsy\n"
             "blue = MetaBrain brain eQTL (this study); purple = deCODE pQTL (prior); "
             "* = relaxed instrument threshold (P<1e-5)", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "panel_B.png"), dpi=300, bbox_inches="tight")
print("panel_B regenerated with", len(fr), "rows")
