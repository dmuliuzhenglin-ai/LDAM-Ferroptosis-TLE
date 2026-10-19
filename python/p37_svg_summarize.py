# -*- coding: utf-8 -*-
"""
Task C (part 2): merge per-sample SVG results, group comparisons, heatmap.

Outputs (results/71_domain_dynamics_svg/):
  svg_morans.csv     gene x sample x Moran's I x permutation P (long format)
  svg_summary.json   per-gene and per-set SE vs CTRL / HS vs CTRL summaries
  figures/v25/Fig_svg_heatmap.png

Run: gsmap_env\\Scripts\\python.exe scripts\\p37_svg_summarize.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results" / "71_domain_dynamics_svg"
PARTS = OUT / "parts"
FIGDIR = BASE / "figures" / "v25"
FIGDIR.mkdir(parents=True, exist_ok=True)

RAT_TP = ["D5", "D10", "D20", "D40"]


def main():
    parts = sorted(PARTS.glob("svg_*.csv"))
    assert len(parts) == 18, f"expected 18 parts, found {len(parts)}"
    df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    df.to_csv(OUT / "svg_morans.csv", index=False)
    print(f"[save] svg_morans.csv rows={len(df)} genes={df['gene'].nunique()}")

    # ---------------------------------------------------------- per-gene summary
    rat = df[df["species"] == "rat"]
    hum = df[df["species"] == "human"]
    per_gene = {}
    for g, sub in df.groupby("gene"):
        r = sub[sub["species"] == "rat"]
        se = r[r["treatment"] == "SE"]
        ct = r[r["treatment"] == "CTRL"]
        h = sub[sub["species"] == "human"]
        hrec = {}
        for _, row in h.iterrows():
            hrec[row["treatment"]] = {"I": round(float(row["morans_i"]), 4),
                                      "p_perm": float(row["p_perm"]),
                                      "sig": bool(row["significant_p005"])}
        per_gene[g] = {
            "gene_sets": sub["gene_sets"].iloc[0],
            "rat_n_detected": int(len(r)),
            "rat_n_sig_SE": int(se["significant_p005"].sum()),
            "rat_n_sig_CTRL": int(ct["significant_p005"].sum()),
            "rat_frac_sig_SE": round(float(se["significant_p005"].mean()), 4)
                               if len(se) else None,
            "rat_frac_sig_CTRL": round(float(ct["significant_p005"].mean()), 4)
                                 if len(ct) else None,
            "rat_mean_I_SE": round(float(se["morans_i"].mean()), 4)
                             if len(se) else None,
            "rat_mean_I_CTRL": round(float(ct["morans_i"].mean()), 4)
                               if len(ct) else None,
            "rat_mean_I": round(float(r["morans_i"].mean()), 4) if len(r) else None,
            "human": hrec,
        }

    # ------------------------------------------------- per-set group comparison
    set_level = {}
    for sname in ["LDAM", "Ferroptosis_core", "Module5", "Hub"]:
        genes = [g for g, v in per_gene.items()
                 if sname in str(v["gene_sets"]).split("|")]
        se_fracs, ct_fracs = [], []
        n_consistent = 0  # sig in >= half of detected rat samples
        for g in genes:
            v = per_gene[g]
            if v["rat_frac_sig_SE"] is not None:
                se_fracs.append(v["rat_frac_sig_SE"])
            if v["rat_frac_sig_CTRL"] is not None:
                ct_fracs.append(v["rat_frac_sig_CTRL"])
            if v["rat_n_detected"] and \
               (v["rat_n_sig_SE"] + v["rat_n_sig_CTRL"]) >= v["rat_n_detected"] / 2:
                n_consistent += 1
        set_level[sname] = {
            "n_genes_tested_rat": len(genes),
            "mean_frac_samples_sig_SE": round(float(np.mean(se_fracs)), 4)
                                        if se_fracs else None,
            "mean_frac_samples_sig_CTRL": round(float(np.mean(ct_fracs)), 4)
                                          if ct_fracs else None,
            "n_genes_sig_in_ge_half_rat_samples": n_consistent,
        }

    # human per-set: fraction sig in HS vs CTRL (n=1 each, descriptive)
    for sname in ["LDAM", "Ferroptosis_core", "Module5", "Hub"]:
        genes = [g for g, v in per_gene.items()
                 if sname in str(v["gene_sets"]).split("|")]
        hs = [per_gene[g]["human"].get("HS", {}).get("sig") for g in genes
              if "HS" in per_gene[g]["human"]]
        ct = [per_gene[g]["human"].get("CTRL", {}).get("sig") for g in genes
              if "CTRL" in per_gene[g]["human"]]
        set_level[sname]["human_frac_sig_HS"] = round(float(np.mean(hs)), 4) if hs else None
        set_level[sname]["human_frac_sig_CTRL"] = round(float(np.mean(ct)), 4) if ct else None
        set_level[sname]["human_n_tested"] = len(hs)

    # per-sample significance counts (traceability)
    per_sample = (df.groupby(["species", "sample", "treatment", "timepoint"])
                    .agg(n_genes=("gene", "count"),
                         n_sig=("significant_p005", "sum"))
                    .reset_index())
    per_sample["frac_sig"] = (per_sample["n_sig"] / per_sample["n_genes"]).round(4)

    summary = {
        "framework": ("Moran's I per gene per sample on log-normalised "
                      "expression (scanpy: min_genes=200, min_cells=3, "
                      "normalize_total 1e4, log1p); KNN k=6 symmetrised "
                      "row-standardised weights; 499 permutations, two-sided "
                      "(|I_perm| >= |I_obs|); seed=42."),
        "n_permutations": 499,
        "alpha": 0.05,
        "rat_group_sizes": {"SE": 8, "CTRL": 8, "per_timepoint_per_group": 2},
        "human_group_sizes": {"HS": 1, "CTRL": 1},
        "per_set": set_level,
        "per_sample": per_sample.to_dict("records"),
        "per_gene": per_gene,
        "limitations": [
            "Spatial transcriptomics expression is strongly spatially "
            "autocorrelated overall, so a large fraction of genes reach "
            "P<0.05; significance here means 'more clustered than random "
            "permutation of spots', not disease association.",
            "Rat: n=2 per group per timepoint (8 SE vs 8 CTRL overall); "
            "human: n=1 HS vs n=1 CTRL. Group comparisons are descriptive.",
            "Permutation P is computed per sample; no cross-sample multiple "
            "testing correction applied to the per-gene calls (FDR noted "
            "where relevant).",
            "Genes not passing QC or lacking a rat ortholog are absent for "
            "that sample ('not detected' is not 'not significant').",
        ],
    }
    with open(OUT / "svg_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("[save] svg_summary.json")
    print(json.dumps(set_level, indent=1))
    make_heatmap(df, per_gene)


def make_heatmap(df, per_gene):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # column order: rat SE D5..D40, rat CTRL D5..D40, human CTRL, human HS
    rat = df[df["species"] == "rat"]
    col_order, col_labels = [], []
    for tr in ["SE", "CTRL"]:
        for tp in RAT_TP:
            for smp in sorted(rat[(rat["treatment"] == tr)
                                  & (rat["timepoint"] == tp)]["sample"].unique()):
                col_order.append(smp)
                col_labels.append(f"{smp.split('_')[0][-4:]}\n{tr} {tp}")
    for smp in ["human_CTRL", "human_HS"]:
        col_order.append(smp)
        col_labels.append(smp.replace("human_", "hum\n"))

    # gene selection: all Hub + Ferroptosis_core detected; then top LDAM and
    # top Module5 by (n_sig_rat desc, mean_I desc), excluding already chosen
    def rank_key(g):
        v = per_gene[g]
        return (v["rat_n_sig_SE"] + v["rat_n_sig_CTRL"], v["rat_mean_I"] or 0)

    def genes_in(sname):
        return [g for g, v in per_gene.items()
                if sname in str(v["gene_sets"]).split("|")
                and v["rat_n_detected"] > 0]

    hub = sorted(genes_in("Hub"), key=rank_key, reverse=True)
    fc = sorted([g for g in genes_in("Ferroptosis_core") if g not in hub],
                key=rank_key, reverse=True)
    ld = sorted([g for g in genes_in("LDAM")
                 if g not in hub and g not in fc], key=rank_key,
                reverse=True)[:14]
    m5 = sorted([g for g in genes_in("Module5")
                 if g not in hub and g not in fc and g not in ld],
                key=rank_key, reverse=True)[:14]
    genes = hub + fc + ld + m5
    print(f"[heatmap] rows={len(genes)} (hub={len(hub)} ferro={len(fc)} "
          f"ldam={len(ld)} m5={len(m5)})")

    piv_I = df.pivot_table(index="gene", columns="sample",
                           values="morans_i", aggfunc="first")
    piv_P = df.pivot_table(index="gene", columns="sample",
                           values="p_perm", aggfunc="first")
    M = piv_I.reindex(index=genes, columns=col_order)
    Pm = piv_P.reindex(index=genes, columns=col_order)

    fig, ax = plt.subplots(figsize=(13.5, 0.32 * len(genes) + 2.2))
    vmax = np.nanmax(np.abs(M.to_numpy()))
    im = ax.imshow(M.to_numpy(), aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_labels, fontsize=6.5, rotation=0)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    # significance stars
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            p = Pm.iat[i, j]
            if pd.notna(p) and p < 0.05:
                ax.text(j, i, "*", ha="center", va="center",
                        fontsize=5, color="black")
    # block separators + labels
    bounds = [0, len(hub), len(hub) + len(fc),
              len(hub) + len(fc) + len(ld), len(genes)]
    names = ["Hub", "Ferroptosis core", "LDAM signature", "Module5"]
    for b in bounds[1:-1]:
        ax.axhline(b - 0.5, color="black", lw=0.8)
    for (b0, b1), nm in zip(zip(bounds[:-1], bounds[1:]), names):
        ax.text(len(col_order) + 0.4, (b0 + b1 - 1) / 2, nm, fontsize=8,
                va="center", fontweight="bold")
    # group separators on x
    ax.axvline(7.5, color="black", lw=0.8)
    ax.axvline(15.5, color="black", lw=0.8)
    ax.text(3.5, -1.3, "Rat SE", ha="center", fontsize=9, fontweight="bold")
    ax.text(11.5, -1.3, "Rat CTRL", ha="center", fontsize=9, fontweight="bold")
    ax.text(16.5, -1.3, "Human", ha="center", fontsize=9, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.08)
    cb.set_label("Moran's I", fontsize=9)
    ax.set_title("Spatial variability (Moran's I) of axis genes across samples\n"
                 "* = permutation P<0.05 (499 perms, KNN k=6), two-sided",
                 fontsize=11, fontweight="bold", pad=26)
    fig.savefig(FIGDIR / "Fig_svg_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("[save] figures/v25/Fig_svg_heatmap.png")


if __name__ == "__main__":
    main()
