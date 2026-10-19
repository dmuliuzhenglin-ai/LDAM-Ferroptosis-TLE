"""Phase 3: figure panels A (MAGMA), B (MR forest), C (LocusCompare LGALS3 cortex)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mr_coloc_utils import genotypes_for, ld_r2

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
FIG = os.path.join(BASE, "figures", "main_figures_v2", "fig9_panels")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.edgecolor": "#333333", "axes.linewidth": 0.8,
                     "figure.facecolor": "white", "axes.facecolor": "white"})

# ---------------- Panel A: MAGMA gene-set bars + hub gene points ----------------
gs = pd.read_csv(os.path.join(BASE, "results", "41_MAGMA", "geneset_competitive_results.csv"))
hub = pd.read_csv(os.path.join(BASE, "results", "41_MAGMA", "hub_genes_genelevel_results.csv"))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.6), gridspec_kw={"width_ratios": [1.1, 1]})
labels = {"LDAM": "LDAM signature", "Ferroptosis": "Ferroptosis", "WGCNA_brown": "WGCNA brown module"}
gs["lab"] = gs["set"].map(labels)
gs = gs.sort_values("P")
colors = ["#2166ac" if q < 0.05 else "#999999" for q in gs["FDR_BH"]]
y = np.arange(len(gs))[::-1]
ax1.barh(y, -np.log10(gs["P"]), color=colors, height=0.6)
for yi, (_, r) in zip(y, gs.iterrows()):
    star = "*" if r["FDR_BH"] < 0.05 else ""
    ax1.text(-np.log10(r["P"]) + 0.03, yi, f"P={r['P']:.3g} q={r['FDR_BH']:.3g}{star}",
             va="center", fontsize=8.5)
ax1.set_yticks(y); ax1.set_yticklabels(gs["lab"])
ax1.axvline(-np.log10(0.05), ls="--", lw=0.8, color="#b2182b")
ax1.text(-np.log10(0.05) + 0.04, -0.42, "P=0.05", color="#b2182b", fontsize=8)
ax1.set_xlabel(r"$-\log_{10}(P)$, competitive gene-set test")
ax1.set_title("MAGMA-style gene-set enrichment\nILAE 2018 all-epilepsy GWAS", fontsize=10)
ax1.set_xlim(0, 3.2)

hub = hub.sort_values("P")
xh = np.arange(len(hub))
ax2.bar(xh, -np.log10(hub["P"]), color="#67a9cf", width=0.6)
for xi, (_, r) in zip(xh, hub.iterrows()):
    ax2.text(xi, -np.log10(r["P"]) + 0.01, f"{r['P']:.2f}", ha="center", fontsize=8)
ax2.set_xticks(xh); ax2.set_xticklabels(hub["symbol"], rotation=30, ha="right")
ax2.set_ylabel(r"$-\log_{10}(P)$")
ax2.set_title("Per-gene MAGMA-style P (hub genes)", fontsize=10)
ax2.set_ylim(0, 1.0)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "panel_A.png"), dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------- Panel B: MR forest ----------------
rows = []
lg = json.load(open(os.path.join(BASE, "results", "42_MR_coloc", "LGALS3_results.json")))
for r in lg:
    mr = r.get("MR:IVW (multiplicative RE)")
    if mr:
        rows.append(dict(label=f"LGALS3 eQTL (MetaBrain cortex, nIV={r['n_instruments']})",
                         b=mr["b"], se=mr["se"], p=mr["p"]))
for f, tag in [("15_MR_LGALS3_EPI_ALL_vs_ieu-b-8.csv", "LGALS3 pQTL (deCODE) vs all epilepsy"),
               ("15_MR_LGALS3_EPI_FOCAL_vs_ieu-b-10.csv", "LGALS3 pQTL (deCODE) vs focal epilepsy"),
               ("15_MR_LGALS3_EPI_HS_vs_ieu-b-14.csv", "LGALS3 pQTL (deCODE) vs focal epilepsy + HS")]:
    d = pd.read_csv(os.path.join(BASE, "results", f))
    d = d[d["method"] == "Inverse variance weighted"].iloc[0]
    rows.append(dict(label=f"{tag} (nIV={int(d['nsnp'])})", b=d["b"], se=d["se"], p=d["pval"]))
fr = pd.DataFrame(rows)
fr["OR"] = np.exp(fr["b"]); fr["lo"] = np.exp(fr["b"] - 1.96 * fr["se"]); fr["hi"] = np.exp(fr["b"] + 1.96 * fr["se"])
fig, ax = plt.subplots(figsize=(7.6, 3.2))
yy = np.arange(len(fr))[::-1]
for yi, (_, r) in zip(yy, fr.iterrows()):
    ax.plot([r["lo"], r["hi"]], [yi, yi], color="#2166ac", lw=1.6, solid_capstyle="round")
    ax.plot(r["OR"], yi, "s", color="#2166ac", ms=7)
    ax.text(1.14, yi, f"OR={r['OR']:.3f} ({r['lo']:.3f}-{r['hi']:.3f}), P={r['p']:.3g}",
            va="center", fontsize=8.5)
ax.axvline(1, ls="--", color="#b2182b", lw=0.9)
ax.set_yticks(yy); ax.set_yticklabels(fr["label"], fontsize=9)
ax.set_xscale("log"); ax.set_xlim(0.88, 1.13)
ax.set_xticks([0.9, 0.95, 1.0, 1.05, 1.1])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel("Odds ratio for epilepsy per unit higher exposure (log scale)")
ax.set_title("Mendelian randomization: LGALS3 (eQTL & pQTL) -> epilepsy", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "panel_B.png"), dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------- Panel C: LocusCompare LGALS3 cortex ----------------
m = pd.read_csv(os.path.join(BASE, "results", "42_MR_coloc", "harmonized_LGALS3_cortex_vs_ILAEall.csv"))
m["lp_eqtl"] = -np.log10(m["p"].clip(1e-300))
m["z_g"] = m["beta_gwas_h"] / m["se_gwas"]
from scipy.stats import norm as _n
m["lp_gwas"] = -np.log10(2 * _n.sf(np.abs(m["z_g"])).clip(1e-300))
lead = m.loc[m["lp_eqtl"].idxmax(), "SNP"]
geno = genotypes_for([lead] + m["SNP"].tolist())
glead = geno[lead]
r2 = np.array([ld_r2(glead, geno[s]) if s in geno else np.nan for s in m["SNP"]])
fig, ax = plt.subplots(figsize=(5.2, 4.6))
sc = ax.scatter(m["lp_gwas"], m["lp_eqtl"], c=r2, cmap="plasma", s=9, vmin=0, vmax=1,
                edgecolors="none", alpha=0.9)
li = m["SNP"] == lead
ax.scatter(m.loc[li, "lp_gwas"], m.loc[li, "lp_eqtl"], facecolors="none",
           edgecolors="#b2182b", s=60, lw=1.4)
ax.annotate(f"lead eQTL SNP {lead}", (m.loc[li, "lp_gwas"].iloc[0], m.loc[li, "lp_eqtl"].iloc[0]),
            textcoords="offset points", xytext=(8, -14), fontsize=8, color="#b2182b")
cb = fig.colorbar(sc, ax=ax, shrink=0.8)
cb.set_label(r"LD $r^2$ with lead SNP (1000G EUR)")
ax.set_xlabel(r"$-\log_{10}(P)$ epilepsy GWAS (ILAE 2018)")
ax.set_ylabel(r"$-\log_{10}(P)$ LGALS3 eQTL (MetaBrain cortex)")
ax.set_title(f"LGALS3 cis region, {len(m)} SNPs\ncoloc PP.H4=0.0084 (no colocalization)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "panel_C.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("lead SNP:", lead, "| panels saved to", FIG)
