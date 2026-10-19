# -*- coding: utf-8 -*-
"""
Task B: Spatial domain composition time dynamics.

Inputs:
  results/67_graphst/domains_per_sample.parquet  (per-spot GraphST domain +
      anatomical annotation for 16 rat samples + 2 human hippocampus samples)
  results/66_spatial_rat/sample_qc.csv           (group / timepoint / spot counts)

Outputs (results/71_domain_dynamics_svg/):
  domain_proportions.csv       per-sample per-domain spot proportions
  domain_timecourse_stats.json SE vs CTRL per timepoint (n=2/group, descriptive)
  figures/v25/Fig_domain_composition_time.png

Run: gsmap_env\\Scripts\\python.exe scripts\\p35_domain_dynamics.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results" / "71_domain_dynamics_svg"
OUT.mkdir(parents=True, exist_ok=True)
FIGDIR = BASE / "figures" / "v25"
FIGDIR.mkdir(parents=True, exist_ok=True)

TP_ORDER = ["D5", "D10", "D20", "D40"]
DOMAIN_ORDER = ["CA1-like", "CA1-like (NETO1+)", "DG-like (PROX1+)",
                "thalamus", "caudoputamen", "white_matter",
                "astro-enriched", "microglia-enriched", "mixed/other"]


def rank_biserial(x, y):
    """x = SE values, y = CTRL values. rbc>0 => SE higher."""
    if len(x) < 2 or len(y) < 2:
        return None, None
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    rbc = 2 * u / (len(x) * len(y)) - 1
    return float(rbc), float(p)


def main():
    df = pd.read_parquet(BASE / "results" / "67_graphst" / "domains_per_sample.parquet")
    qc = pd.read_csv(BASE / "results" / "66_spatial_rat" / "sample_qc.csv")

    # ------------------------------------------------ per-sample proportions
    counts = (df.groupby(["sample", "graphst_annotation"]).size()
                .rename("n").reset_index())
    totals = df.groupby("sample").size().rename("n_total").reset_index()
    meta = (df.groupby("sample")[["species", "treatment", "timepoint"]]
              .first().reset_index())
    prop = counts.merge(totals, on="sample").merge(meta, on="sample")
    prop["proportion"] = prop["n"] / prop["n_total"]
    prop = prop.rename(columns={"graphst_annotation": "domain"})
    # complete sample x domain grid: absent domain in a sample = 0 spots
    full_idx = pd.MultiIndex.from_product(
        [meta["sample"], sorted(df["graphst_annotation"].unique())],
        names=["sample", "domain"])
    prop = (prop.set_index(["sample", "domain"]).reindex(full_idx)
                .reset_index())
    prop = prop.merge(meta, on="sample", suffixes=("", "_m"))
    for c in ["species", "treatment", "timepoint"]:
        prop[c] = prop[c].fillna(prop[c + "_m"])
    prop = prop.merge(totals, on="sample", suffixes=("", "_t"))
    prop["n"] = prop["n"].fillna(0).astype(int)
    prop["n_total"] = prop["n_total"].fillna(prop["n_total_t"]).astype(int)
    prop["proportion"] = prop["proportion"].fillna(0.0)
    prop = prop[["species", "sample", "treatment", "timepoint", "domain",
                 "n", "n_total", "proportion"]]
    prop.to_csv(OUT / "domain_proportions.csv", index=False)
    print(f"[save] domain_proportions.csv rows={len(prop)}")

    rat = prop[prop["species"] == "rat"]
    hum = prop[prop["species"] == "human"]
    domains = [d for d in DOMAIN_ORDER if d in set(rat["domain"])]
    for d in sorted(set(rat["domain"]) - set(domains)):
        domains.append(d)

    # -------------------------------------- rat SE vs CTRL x timepoint stats
    rat_stats = {}
    for dom in domains:
        dom_rec = {}
        for tp in TP_ORDER:
            se = rat[(rat["domain"] == dom) & (rat["timepoint"] == tp)
                     & (rat["treatment"] == "SE")]["proportion"].tolist()
            ct = rat[(rat["domain"] == dom) & (rat["timepoint"] == tp)
                     & (rat["treatment"] == "CTRL")]["proportion"].tolist()
            rbc, p = rank_biserial(se, ct)
            dom_rec[tp] = {
                "SE_props": [round(v, 5) for v in se],
                "CTRL_props": [round(v, 5) for v in ct],
                "mean_SE": float(np.mean(se)) if se else None,
                "mean_CTRL": float(np.mean(ct)) if ct else None,
                "delta_SE_minus_CTRL": (float(np.mean(se) - np.mean(ct))
                                        if se and ct else None),
                "rank_biserial_r": rbc,
                "mannwhitney_p_descriptive": p,
                "note": "n=2 per group; descriptive only",
            }
        rat_stats[dom] = dom_rec

    # per-treatment timecourse direction (mean across the 2 replicates)
    timecourse_dir = {}
    for dom in domains:
        rec = {}
        for tr in ["SE", "CTRL"]:
            means = []
            for tp in TP_ORDER:
                v = rat_stats[dom][tp][f"mean_{tr}"]
                means.append(v if v is not None else np.nan)
            means = np.array(means, float)
            rec[tr] = {"means_by_timepoint": dict(zip(TP_ORDER,
                           [None if np.isnan(m) else round(float(m), 5)
                            for m in means])),
                       "delta_D40_minus_D5": (round(float(means[-1] - means[0]), 5)
                                              if not np.isnan(means).any() else None)}
        timecourse_dir[dom] = rec

    # human descriptive
    hum_rec = {}
    for smp, sub in hum.groupby("sample"):
        hum_rec[smp] = {r["domain"]: round(float(r["proportion"]), 5)
                        for _, r in sub.iterrows()}
        hum_rec[smp]["_n_total"] = int(sub["n_total"].iloc[0])

    qc_small = qc[["sample", "treatment", "timepoint", "n_spots_in_tissue",
                   "n_spots_after_qc"]].to_dict("records")
    stats_json = {
        "framework": ("Per-sample spot proportions of GraphST-annotated "
                      "domains; proportions used (not absolute counts) to "
                      "mitigate differing tissue coverage / capture size."),
        "n_per_group_per_timepoint": 2,
        "rat_se_vs_ctrl_by_timepoint": rat_stats,
        "rat_timecourse_means": timecourse_dir,
        "human_descriptive": hum_rec,
        "rat_sample_spot_counts": qc_small,
        "limitations": [
            "n=2 biological replicates per group per timepoint: all "
            "between-group contrasts are descriptive (effect sizes / "
            "exploratory Mann-Whitney only, no inferential power).",
            "Sections differ in total captured spots and anatomical "
            "coverage; proportions (compositional data) are used, so an "
            "apparent expansion of one domain can reflect shrinkage of "
            "another (compositionality confound).",
            "GraphST domains are data-driven clusters annotated by marker "
            "genes; anatomical labels are approximate.",
            "Human comparison: 1 HS vs 1 CTRL (GSE288888), descriptive only.",
        ],
    }
    with open(OUT / "domain_timecourse_stats.json", "w") as f:
        json.dump(stats_json, f, indent=2)
    print("[save] domain_timecourse_stats.json")

    make_figure(rat, hum, domains)


def make_figure(rat, hum, domains):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    cmap = {"CA1-like": "#4C72B0", "CA1-like (NETO1+)": "#8DA9DB",
            "DG-like (PROX1+)": "#55A868", "thalamus": "#C44E52",
            "caudoputamen": "#DD8452", "white_matter": "#8172B3",
            "astro-enriched": "#CCB974", "microglia-enriched": "#D62728",
            "mixed/other": "#BBBBBB"}

    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.35, 1, 1],
                          hspace=0.55, wspace=0.38,
                          left=0.07, right=0.90, top=0.93, bottom=0.06)

    # ---------------- Panel A: stacked bars of mean proportions
    axA = fig.add_subplot(gs[0, :2])
    combos = [(tr, tp) for tp in TP_ORDER for tr in ["SE", "CTRL"]]
    x = np.arange(len(combos))
    bottom = np.zeros(len(combos))
    for dom in domains:
        vals = []
        for tr, tp in combos:
            sub = rat[(rat["domain"] == dom) & (rat["treatment"] == tr)
                      & (rat["timepoint"] == tp)]["proportion"]
            vals.append(sub.mean() if len(sub) else 0.0)
        vals = np.array(vals)
        axA.bar(x, vals, bottom=bottom, width=0.75,
                color=cmap.get(dom, "#999999"), label=dom,
                edgecolor="white", linewidth=0.3)
        bottom += vals
    axA.set_xticks(x)
    axA.set_xticklabels([f"{tr}\n{tp}" for tr, tp in combos], fontsize=8)
    axA.set_ylabel("Mean spot proportion")
    axA.set_title("A  Rat domain composition (mean of n=2 replicates)",
                  fontsize=11, fontweight="bold", loc="left")
    axA.set_ylim(0, 1.0)
    axA.spines[["top", "right"]].set_visible(False)

    # ---------------- Panel B: human stacked bars (descriptive)
    axB = fig.add_subplot(gs[0, 2])
    hum_piv = (hum.pivot_table(index="sample", columns="domain",
                               values="proportion", aggfunc="sum")
                 .reindex(columns=[d for d in domains]))
    hum_piv = hum_piv.fillna(0.0)
    xb = np.arange(len(hum_piv))
    bot = np.zeros(len(hum_piv))
    for dom in domains:
        if dom in hum_piv.columns:
            v = hum_piv[dom].to_numpy()
            axB.bar(xb, v, bottom=bot, width=0.6, color=cmap.get(dom, "#999999"),
                    edgecolor="white", linewidth=0.3)
            bot += v
    axB.set_xticks(xb)
    axB.set_xticklabels(hum_piv.index, fontsize=8)
    axB.set_ylim(0, 1.0)
    axB.set_title("B  Human hippocampus\n(GSE288888, n=1 each)",
                  fontsize=11, fontweight="bold", loc="left")
    axB.spines[["top", "right"]].set_visible(False)

    # legend for A/B
    axL = fig.add_subplot(gs[0, 3])
    axL.axis("off")
    handles = [Patch(facecolor=cmap.get(d, "#999999"), label=d) for d in domains]
    axL.legend(handles=handles, loc="center left", fontsize=9,
               frameon=False, title="Domain", title_fontsize=10)

    # ---------------- Panels C: per-domain dot timecourse (3x? grid)
    focus = [d for d in ["microglia-enriched", "astro-enriched", "CA1-like",
                         "CA1-like (NETO1+)", "thalamus", "caudoputamen",
                         "white_matter", "mixed/other"] if d in domains]
    letters = "CDEFGHIJ"
    col_map = {"SE": "#D62728", "CTRL": "#1F77B4"}
    for i, dom in enumerate(focus):
        r, c = divmod(i, 4)
        ax = fig.add_subplot(gs[r + 1, c])
        xticks = np.arange(len(TP_ORDER))
        for tr in ["SE", "CTRL"]:
            for tp_i, tp in enumerate(TP_ORDER):
                sub = rat[(rat["domain"] == dom) & (rat["treatment"] == tr)
                          & (rat["timepoint"] == tp)]["proportion"]
                jitter = -0.13 if tr == "SE" else 0.13
                ax.scatter(np.full(len(sub), tp_i + jitter), sub,
                           color=col_map[tr], s=22, alpha=0.9, zorder=3,
                           edgecolor="white", linewidth=0.4)
            means = [rat[(rat["domain"] == dom) & (rat["treatment"] == tr)
                         & (rat["timepoint"] == tp)]["proportion"].mean()
                     for tp in TP_ORDER]
            off = -0.13 if tr == "SE" else 0.13
            ax.plot(xticks + off, means, color=col_map[tr], lw=1.6,
                    marker="o", ms=4, label=tr, zorder=2)
        ax.set_xticks(xticks)
        ax.set_xticklabels(TP_ORDER, fontsize=8)
        ax.set_title(f"{letters[i]}  {dom}", fontsize=10, fontweight="bold",
                     loc="left")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
        if c == 0:
            ax.set_ylabel("Spot proportion", fontsize=9)
        if i == 0:
            ax.legend(fontsize=8, frameon=False, loc="best")
    fig.suptitle("Spatial domain composition dynamics — rat SE vs CTRL "
                 "timecourse (n=2/group) and human hippocampus",
                 fontsize=13, fontweight="bold")
    out = FIGDIR / "Fig_domain_composition_time.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[save] {out}")


if __name__ == "__main__":
    main()
