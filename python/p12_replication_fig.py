"""Replication figure: ILAE2018 vs ILAE2023 across MAGMA geneset, scDRS, gsMap + fine-mapping summary."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
FIGD = os.path.join(BASE, "figures", "replication_v24")
os.makedirs(FIGD, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 9})

# ---------- Panel A: MAGMA competitive enrichment ----------
rows = []
src = [("ILAE2018 all epilepsy", os.path.join(BASE, "results", "41_MAGMA", "geneset_competitive_results.csv"))]
for ph, lab in [("all_epilepsy", "ILAE2023 all epilepsy"), ("focal_epilepsy", "ILAE2023 focal"),
                ("focal_HS", "ILAE2023 focal-HS")]:
    src.append((lab, os.path.join(BASE, "results", "51_ILAE2023_MAGMA", ph, "geneset_competitive_results.csv")))
for lab, f in src:
    d = pd.read_csv(f)
    for _, r in d.iterrows():
        if r["set"] in ("LDAM", "Ferroptosis"):
            rows.append(dict(dataset=lab, geneset=r["set"], beta=r["BETA"], se=r["SE"], p=r["P"], fdr=r["FDR_BH"]))
A = pd.DataFrame(rows)
A.to_csv(os.path.join(FIGD, "panelA_magma_replication.csv"), index=False)

# ---------- Panel B: scDRS cell-type t-stats ----------
b18 = pd.read_csv(os.path.join(BASE, "results", "43_scDRS", "downstream_celltype_association.csv"))
b18["dataset"] = "ILAE2018 all"
bs = [b18]
for t, lab in [("ILAE3_allEpilepsy", "ILAE2023 all"), ("ILAE3_focalEpilepsy", "ILAE2023 focal"),
               ("ILAE3_focalHS", "ILAE2023 focal-HS")]:
    d = pd.read_csv(os.path.join(BASE, "results", "52_ILAE2023_scDRS", f"downstream_celltype_association_{t}.csv"))
    d["dataset"] = lab
    bs.append(d)
B = pd.concat(bs)
Bp = B.pivot(index="celltype", columns="dataset", values="tstat")[["ILAE2018 all", "ILAE2023 all", "ILAE2023 focal", "ILAE2023 focal-HS"]]
Bp = Bp.loc[["Microglia", "Astrocyte", "Oligodendrocyte", "OPC", "ExcNeuron", "InhNeuron", "Endothelial"]]
Bp.to_csv(os.path.join(FIGD, "panelB_scdrs_tstat.csv"))

# ---------- Panel C: gsMap Cauchy ----------
crows = []
f18 = os.path.join(BASE, "results", "50_gsMap", "quick_mode_Epilepsy",
                   "V1_Mouse_Brain_Sagittal_Posterior_Section", "cauchy_combination",
                   "V1_Mouse_Brain_Sagittal_Posterior_Section_Epilepsy.Cauchy.csv.gz")
if os.path.exists(f18):
    d = pd.read_csv(f18); d["dataset"] = "ILAE2018 all"; crows.append(d)
for t, lab in [("all_epilepsy", "ILAE2023 all"), ("focal_epilepsy", "ILAE2023 focal"), ("focal_HS", "ILAE2023 focal-HS")]:
    f = os.path.join(BASE, "results", "53_ILAE2023_gsMap", f"quick_mode_{t}",
                     "V1_Mouse_Brain_Sagittal_Posterior_Section", "cauchy_combination",
                     f"V1_Mouse_Brain_Sagittal_Posterior_Section_{t}.Cauchy.csv.gz")
    d = pd.read_csv(f); d["dataset"] = lab; crows.append(d)
C = pd.concat(crows)
C["mlog10p"] = -np.log10(C["p_cauchy"])
Cp = C.pivot(index="annotation", columns="dataset", values="mlog10p")
DOM = {1: "D1 white matter", 2: "D2 cortex L2/3-6", 3: "D3 cortex L5/6", 4: "D4 hippocampus", 5: "D5 dentate gyrus", 6: "D6 white matter"}
Cp.index = [f"D{i} ({DOM.get(i,'')})" for i in Cp.index]
Cp.to_csv(os.path.join(FIGD, "panelC_gsmap_cauchy.csv"))

# ---------- Panel D: fine-mapping ----------
F = pd.read_csv(os.path.join(BASE, "results", "55_finemap", "finemap_summary_ILAE2018_all.csv"))
Fd = F[F["lead_p"] < 5e-6].head(12)

# ---------- draw ----------
fig, axes = plt.subplots(1, 4, figsize=(15, 4.2), gridspec_kw={"width_ratios": [1.1, 1.3, 1.1, 1.3]})

ax = axes[0]
ys, xs, xlo, xhi, cs, labs = [], [], [], [], [], []
y = 0
for gs in ["LDAM", "Ferroptosis"]:
    for lab in ["ILAE2018 all epilepsy", "ILAE2023 all epilepsy", "ILAE2023 focal", "ILAE2023 focal-HS"]:
        r = A[(A["geneset"] == gs) & (A["dataset"] == lab)]
        if len(r) == 0: continue
        r = r.iloc[0]
        ys.append(y); xs.append(r["beta"]); xlo.append(r["beta"] - 1.96 * r["se"]); xhi.append(r["beta"] + 1.96 * r["se"])
        cs.append("#c0392b" if r["fdr"] < 0.05 else ("#e67e22" if r["p"] < 0.05 else "#7f8c8d"))
        star = "*" if r["fdr"] < 0.05 else ("+" if r["p"] < 0.1 else "")
        labs.append(f"{lab.replace(' epilepsy','')}{star}")
        y += 1
    y += 0.8
ax.errorbar(xs, ys, xerr=[np.array(xs) - np.array(xlo), np.array(xhi) - np.array(xs)], fmt="none", ecolor="#bbbbbb", zorder=1)
ax.scatter(xs, ys, c=cs, s=45, zorder=2)
ax.axvline(0, color="k", lw=0.6, ls="--")
ax.set_yticks(ys); ax.set_yticklabels(labs, fontsize=7.5)
ax.set_xlabel("MAGMA competitive enrichment (beta)")
ax.set_title("A  Gene-set enrichment replication", loc="left", fontweight="bold")
# group separators
ax.text(ax.get_xlim()[1], 1.4, "LDAM", fontsize=8, fontweight="bold", ha="right")
ax.text(ax.get_xlim()[1], 6.2, "Ferroptosis", fontsize=8, fontweight="bold", ha="right")

ax = axes[1]
im = ax.imshow(Bp.values, cmap="RdBu_r", vmin=-40, vmax=40, aspect="auto")
ax.set_xticks(range(len(Bp.columns))); ax.set_xticklabels(Bp.columns, rotation=30, ha="right", fontsize=7.5)
ax.set_yticks(range(len(Bp.index))); ax.set_yticklabels(Bp.index, fontsize=7.5)
for i in range(Bp.shape[0]):
    for j in range(Bp.shape[1]):
        v = Bp.values[i, j]
        ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.5,
                color="white" if abs(v) > 22 else "black")
ax.set_title("B  scDRS cell-type association (t)", loc="left", fontweight="bold")
fig.colorbar(im, ax=ax, shrink=0.6)

ax = axes[2]
x = np.arange(len(Cp.index)); w = 0.2
for k, col in enumerate(Cp.columns):
    ax.bar(x + (k - 1.5) * w, Cp[col], width=w, label=col, alpha=0.85)
ax.axhline(-np.log10(0.05), color="k", ls="--", lw=0.7)
ax.set_xticks(x); ax.set_xticklabels(Cp.index, rotation=30, ha="right", fontsize=7)
ax.set_ylabel("-log10 Cauchy P")
ax.legend(fontsize=6, frameon=False)
ax.set_title("C  gsMap spatial enrichment", loc="left", fontweight="bold")

ax = axes[3]
Fd = Fd.sort_values("lead_p")
ax.scatter(-np.log10(Fd["lead_p"]), range(len(Fd)), s=20 + 200 * Fd["lead_pp"], c="#2c3e50", alpha=0.75)
for i, (_, r) in enumerate(Fd.iterrows()):
    ax.text(-np.log10(r["lead_p"]) + 0.05, i, f'{r["gene"]} (PP={r["lead_pp"]:.2f}, CS={r["cs95_size"]})', fontsize=6.5, va="center")
ax.set_yticks([])
ax.set_xlabel("lead SNP -log10 P")
ax.set_xlim(left=min(-np.log10(Fd["lead_p"])) - 0.2)
ax.set_title("D  ABF fine-mapping (top loci)", loc="left", fontweight="bold")

plt.tight_layout()
fig.savefig(os.path.join(FIGD, "Figure_replication_v24.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIGD, "Figure_replication_v24.pdf"), bbox_inches="tight")
print("figure saved")
