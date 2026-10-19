"""p22: TWAS figures (managed python).
Fig_TWAS_manhattan.png  - hippocampus x focal-HS Manhattan, FDR<0.05 labeled
                          (top-5 nominal if none pass FDR)
Fig_TWAS_overlap.png    - overlap of TWAS genes with LDAM+Ferroptosis axis
                          gene set, Fisher exact test, OR and P annotated
"""
import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as st

BASE = os.getcwd()
OUT = os.path.join(BASE, "results", "61_twas")
FIG = os.path.join(BASE, "figures", "v25")
os.makedirs(FIG, exist_ok=True)
GTF = os.path.join(BASE, "results", "50_gsMap", "gsMap_resource",
                   "genome_annotation", "gtf",
                   "gencode.v46lift37.basic.annotation.gtf")

def ens_maps():
    e2s, s2e = {}, {}
    with open(GTF, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            a = c[8]
            gid = a.split('gene_id "')[1].split('"')[0].split(".")[0].split("_")[0]
            gn = a.split('gene_name "')[1].split('"')[0] if 'gene_name "' in a else gid
            e2s[gid] = gn
            s2e.setdefault(gn, gid)
    return e2s, s2e

e2s, s2e = ens_maps()

# ---------------- Manhattan: Hippocampus x focal-HS ----------------
R = pd.read_csv(os.path.join(OUT, "TWAS_Hippocampus_fhs.csv"))
R["sym"] = [e2s.get(g.split(".")[0], g) for g in R.gene]
R["mid"] = (R.p0 + R.p1) / 2
chr_lens = R.groupby("chr").mid.max()
offsets = {}
cum = 0
for c in sorted(chr_lens.index):
    offsets[c] = cum
    cum += chr_lens[c] + 5e6
R["x"] = [m + offsets[c] for m, c in zip(R.mid, R.chr)]
xticks = [offsets[c] + chr_lens[c] / 2 for c in sorted(chr_lens.index)]

fig, ax = plt.subplots(figsize=(13, 4.8))
colors = ["#4C72B0", "#55A868"]
for i, (c, grp) in enumerate(R.groupby("chr")):
    ax.scatter(grp.x, -np.log10(grp.P), s=9, color=colors[i % 2],
               alpha=0.75, linewidths=0)
bonf = 0.05 / len(R)
ax.axhline(-np.log10(bonf), color="#C44E52", ls="--", lw=1,
           label=f"Bonferroni 0.05/{len(R)} (P={bonf:.1e})")
sig = R[R.FDR < 0.05]
if len(sig):
    lab = sig
    note = f"FDR<0.05: {len(sig)} genes"
else:
    lab = R.nsmallest(5, "P")
    note = "no gene passes FDR<0.05; top 5 nominal labelled"
for _, r in lab.iterrows():
    ax.annotate(r.sym, (r.x, -np.log10(r.P)), textcoords="offset points",
                xytext=(0, 6), ha="center", fontsize=8.5, fontweight="bold")
    ax.scatter([r.x], [-np.log10(r.P)], s=22, color="#C44E52", zorder=5)
ax.set_xticks(xticks)
ax.set_xticklabels([str(c) for c in sorted(chr_lens.index)], fontsize=8)
ax.set_xlabel("Chromosome")
ax.set_ylabel(r"$-\log_{10}(P)$")
ax.set_title("TWAS Manhattan - GTEx v8 Brain Hippocampus × ILAE 2023 focal epilepsy with HS\n"
             + note, fontsize=11)
ax.legend(frameon=False, fontsize=8, loc="upper right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "Fig_TWAS_manhattan.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("manhattan done;", note)

# ---------------- overlap with LDAM+Ferroptosis axis genes ----------------
axis_syms = set(pd.read_csv(os.path.join(BASE, "genesets", "LDAM_Ferroptosis_union.txt"),
                            header=None)[0].str.strip())
if len(axis_syms) == 0:  # fallback: cismr gene column
    axis_syms = set(pd.read_csv(os.path.join(BASE, "results", "60_axis_cismr",
                                             "axis_cismr_results.csv")).gene.unique())
axis_ensg = {s2e[s] for s in axis_syms if s in s2e}
print("axis genes:", len(axis_syms), "mapped to ENSG:", len(axis_ensg))

R["ensg"] = [g.split(".")[0] for g in R.gene]
tested = set(R.ensg)
axis_t = tested & axis_ensg

panels = []
for label, sigset in [("FDR < 0.05", set(R.loc[R.FDR < 0.05, "ensg"])),
                      ("nominal P < 0.05", set(R.loc[R.P < 0.05, "ensg"]))]:
    a = len(sigset & axis_t)
    b = len(axis_t - sigset)
    c = len(sigset - axis_t)
    d = len(tested - sigset - axis_t)
    if a + c == 0:
        OR, P = np.nan, np.nan
    else:
        OR, P = st.fisher_exact([[a, b], [c, d]])
    panels.append((label, a, b, c, d, OR, P))
    print(label, "table:", [[a, b], [c, d]], "OR:", OR, "P:", P)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
for ax, (label, a, b, c, d, OR, P) in zip(axes, panels):
    cats = ["Axis genes\nTWAS+", "Axis genes\nTWAS-", "Other genes\nTWAS+", "Other genes\nTWAS-"]
    vals = [a, b, c, d]
    bars = ax.bar(cats, vals, color=["#C44E52", "#DD8452", "#4C72B0", "#8C8C8C"], width=0.65)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v, str(v), ha="center",
                va="bottom", fontsize=10, fontweight="bold")
    or_txt = f"OR = {OR:.2f}, Fisher P = {P:.3f}" if not np.isnan(OR) else "OR undefined (0 TWAS+ genes)"
    ax.set_title(f"TWAS+ = {label}\n{or_txt}", fontsize=10)
    ax.set_ylabel("gene count")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=8)
fig.suptitle("Overlap of hippocampus × focal-HS TWAS genes with LDAM+Ferroptosis axis genes\n"
             f"(universe: {len(tested)} tested genes; axis genes tested: {len(axis_t)})",
             fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "Fig_TWAS_overlap.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("overlap fig done")
