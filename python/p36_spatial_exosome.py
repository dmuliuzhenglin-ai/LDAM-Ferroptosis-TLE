# -*- coding: utf-8 -*-
"""
p36_spatial_exosome.py
Spatial co-localisation test: exosome secretion/uptake machinery x COMMOT Gal-3
broadcast field.

Hypothesis: the tissue-level broadcast of signal-peptide-less Gal-3 (LGALS3,
a validated exosome cargo; results/40_exosome_cargo) is mediated by
extracellular vesicles. Testable prediction: the spatial distribution of
exosome secretion / uptake machinery scores should co-localise with the
COMMOT Gal-3 sender/receiver fields and be temporally synchronised.

Gene-set provenance (traceable): the 7 exosome modules + 2 reference sets are
copied verbatim from R/27_exosome_function.R lines 39-54 (the same definitions
used by analyses 27 / 28 / 36).

  Human: GSE288888, GSM8778056 (CTRL) / GSM8778057 (MTLE-HS), n=1 vs n=1.
  Rat:   GSE301260, Li-pilocarpine SE, D5/D10/D20/D40 x SE/CTRL x 2 = 16.

Statistical discipline:
  - human n=1/1 and rat n=2/group -> effect sizes, spot-level tests descriptive;
  - the PRIMARY inferential test is sample-level Spearman across the 16 rat
    samples (module mean vs per-sample Gal-3 receiver range/intensity).

Run:  gsmap_env\\Scripts\\python.exe scripts\\p36_spatial_exosome.py [--figs-only]
Outputs: results/72_spatial_exosome/ + figures/v25/Fig_spatial_exosome.png
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[1]
DATA_H = ROOT / "data" / "spatial" / "GSE288888"
DATA_R = ROOT / "data" / "spatial" / "GSE301260"
OUT = ROOT / "results" / "72_spatial_exosome"
CACHE = OUT / "score_cache"
FIGD = ROOT / "figures" / "v25"
for d in (OUT, CACHE, FIGD):
    d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Gene sets: verbatim copy of R/27_exosome_function.R lines 39-54.
# --------------------------------------------------------------------------
GENE_SETS = {
    "ESCRT_core": ["TSG101", "VPS4A", "VPS4B", "VPS25", "VPS36", "HGS", "STAM",
                   "STAM2", "CHMP1A", "CHMP1B", "CHMP2A", "CHMP2B", "CHMP3",
                   "CHMP4A", "CHMP4B", "CHMP4C", "CHMP5", "CHMP6", "CHMP7",
                   "PDCD6IP"],
    "Rab_secretory": ["RAB27A", "RAB27B", "RAB11A", "RAB11B", "RAB35", "RAB8A",
                      "RAB8B", "RAB7A"],
    "Ceramide_nSMase": ["SMPD2", "SMPD3", "SMPD4"],
    "Cargo_loading": ["YBX1", "SYNCRIP", "HNRNPA2B1", "FMR1", "LIN28B",
                      "IGF2BP1"],
    "SNARE_fusion": ["SNAP23", "STX4", "STX6", "STX7", "VAMP3", "VAMP7",
                     "VTI1B"],
    "Syntenin_ectosome": ["SDCBP", "SDCBP2", "PDCD6IP", "ARRDC1", "VPS4A",
                          "VPS4B"],
    "Exosome_uptake": ["CD9", "CD81", "ICAM1", "ITGB1", "ITGAV", "SDC1",
                       "SDC4", "HSPG2"],
    "LDAM": ["PLIN2", "DGAT1", "DGAT2", "ACAT1", "ACSL1", "CD36", "FABP3",
             "HILPDA", "SOAT1", "LPIN1"],
    "Ferroptosis": ["GPX4", "ACSL4", "LPCAT3", "ALOX15", "AIFM2", "TFRC",
                    "FTH1", "FTL", "SLC7A11", "NCOA4"],
}
MODULES = ["ESCRT_core", "Rab_secretory", "Ceramide_nSMase", "Cargo_loading",
           "SNARE_fusion", "Syntenin_ectosome", "Exosome_uptake"]
REFS = ["LDAM", "Ferroptosis"]

HUMAN_SAMPLES = {"CTRL": "GSM8778056_A18-43", "HS": "GSM8778057_276A"}


# ---------------------------------------------------------------- loading
def read_positions(path):
    with gzip.open(path, "rt") as f:
        first = f.readline()
    cols = ["barcode", "in_tissue", "array_row", "array_col",
            "pxl_row_in_fullres", "pxl_col_in_fullres"]
    if first.lower().startswith("barcode"):
        pos = pd.read_csv(path)
        pos.columns = [c.strip() for c in pos.columns]
    else:
        pos = pd.read_csv(path, header=None, names=cols)
    return pos


def load_sample(prefix, pos_suffix):
    with gzip.open(f"{prefix}_matrix.mtx.gz", "rb") as f:
        mat = scipy.io.mmread(f)
    X = sparse.csr_matrix(mat).T.tocsr()
    feats = pd.read_csv(f"{prefix}_features.tsv.gz", sep="\t", header=None,
                        names=["id", "symbol", "ftype"])
    bc = pd.read_csv(f"{prefix}_barcodes.tsv.gz", header=None,
                     names=["barcode"])
    pos = read_positions(f"{prefix}_{pos_suffix}.gz")
    ad = sc.AnnData(X=X, obs=bc.set_index("barcode"),
                    var=feats.set_index("symbol"))
    ad.var_names_make_unique()
    pos = pos.set_index("barcode")
    ad.obs = ad.obs.join(pos[["in_tissue", "array_row", "array_col",
                              "pxl_row_in_fullres", "pxl_col_in_fullres"]])
    ad = ad[ad.obs["in_tissue"] == 1].copy()
    ad.obsm["spatial"] = ad.obs[["pxl_col_in_fullres",
                                 "pxl_row_in_fullres"]].to_numpy(float)
    return ad


def qc_normalize(ad, mt_regex):
    ad.var["mt"] = ad.var_names.str.match(mt_regex)
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    return ad


# ------------------------------------------------------- rat ortholog map
def rat_map():
    """Human symbol -> rat symbol: cached mygene/title-case map first,
    title-case fallback for genes not covered (provenance recorded)."""
    om = json.load(open(ROOT / "results" / "66_spatial_rat" /
                        "ortholog_map.json"))
    prov, out = {}, {}
    all_genes = sorted({g for gs in GENE_SETS.values() for g in gs})
    for g in all_genes:
        if g in om and om[g]:
            out[g], prov[g] = om[g], "ortholog_map.json"
        else:
            # first letter upper, rest lower (TSG101 -> Tsg101), matching
            # the to_mouse() convention in R/27_exosome_function.R line 57
            out[g], prov[g] = g[0].upper() + g[1:].lower(), \
                "first_upper_rest_lower_fallback"
    return out, prov


# ---------------------------------------------------------------- scoring
def score_one(ad, sets, prefix="score_"):
    avail = set(ad.var_names)
    det = {}
    for name, gs in sets.items():
        gs_use = [g for g in gs if g in avail]
        det[name] = {"n_input": len(gs), "n_detected": len(gs_use),
                     "missing": [g for g in gs if g not in avail]}
        if len(gs_use) >= 2:
            sc.tl.score_genes(ad, gs_use, ctrl_size=50,
                              score_name=f"{prefix}{name}", random_state=0)
        else:
            ad.obs[f"{prefix}{name}"] = np.nan
    return det


def rank_biserial(x, y):
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    return float(u), float(p), float(2 * u / (len(x) * len(y)) - 1)


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


# ================================================================== compute
def compute():
    rat_orth, rat_prov = rat_map()
    rat_sets = {k: [rat_orth[g] for g in v] for k, v in GENE_SETS.items()}
    json.dump(rat_prov, open(OUT / "rat_ortholog_provenance.json", "w"),
              indent=2)

    commot = pd.read_parquet(ROOT / "results" / "70_commot_spatial" /
                             "lr_signal_per_spot.parquet")
    commot = commot[["sample", "barcode", "gal3_sender_total",
                     "gal3_receiver_total"]].rename(columns={"sample": "gsm"})

    # ---------------- human: merged common-reference scoring (as p30) ----
    cache_h = CACHE / "human_scores.parquet"
    if cache_h.exists():
        hdf = pd.read_parquet(cache_h).reset_index(drop=True)
        det_h = json.load(open(CACHE / "human_detected.json"))
    else:
        ads = {t: qc_normalize(load_sample(DATA_H / s, "tissue_positions.csv"),
                               r"(?i)^MT-")
               for t, s in HUMAN_SAMPLES.items()}
        for t, ad in ads.items():
            ad.obs["sample"] = t
            print(f"[load] human {t}: {ad.n_obs} spots")
        comb = sc.concat(ads, label="batch", join="inner", index_unique="-")
        det_h = score_one(comb, GENE_SETS)
        cols = (["sample", "array_row", "array_col", "pxl_row_in_fullres",
                 "pxl_col_in_fullres"]
                + [f"score_{m}" for m in MODULES + REFS])
        hdf = comb.obs[cols].copy()
        hdf = hdf.reset_index(drop=True)
        hdf.insert(0, "barcode",
                   [b.rsplit("-", 1)[0] for b in comb.obs_names])
        hdf["gsm"] = hdf["sample"].map(
            dict(zip(HUMAN_SAMPLES.keys(), HUMAN_SAMPLES.values())))
        hdf.to_parquet(cache_h)
        json.dump(det_h, open(CACHE / "human_detected.json", "w"), indent=2)
    n_h = hdf.merge(commot, on=["gsm", "barcode"], how="inner")
    print(f"[join] human spots matched to COMMOT: "
          f"{len(n_h)}/{len(hdf)}")
    n_h.to_parquet(OUT / "human_spot_table.parquet")

    # ---------------- rat: per-sample scoring (as p31) -------------------
    qc = pd.read_csv(ROOT / "results" / "66_spatial_rat" / "sample_qc.csv")
    frames, det_r = [], {}
    for _, row in qc.iterrows():
        smp = row["sample"]
        fc = CACHE / f"rat_{smp}.parquet"
        if fc.exists():
            df = pd.read_parquet(fc)
            det_r[smp] = json.load(open(CACHE / f"rat_{smp}_detected.json"))
        else:
            ad = qc_normalize(load_sample(DATA_R / smp,
                                          "tissue_positions_list.csv"),
                              r"(?i)^mt[-_]")
            print(f"[load] rat {smp}: {ad.n_obs} spots")
            det_r[smp] = score_one(ad, rat_sets)
            cols = (["array_row", "array_col", "pxl_row_in_fullres",
                     "pxl_col_in_fullres"]
                    + [f"score_{m}" for m in MODULES + REFS])
            df = ad.obs[cols].copy()
            df = df.reset_index(drop=True)
            df.insert(0, "barcode", ad.obs_names)
            df["sample"] = smp
            df["treatment"] = row["treatment"]
            df["timepoint"] = row["timepoint"]
            df.to_parquet(fc)
            json.dump(det_r[smp],
                      open(CACHE / f"rat_{smp}_detected.json", "w"), indent=2)
        frames.append(df)
    rdf = pd.concat(frames, ignore_index=True)
    rdf["gsm"] = rdf["sample"]
    n_r = rdf.merge(commot, on=["gsm", "barcode"], how="inner")
    print(f"[join] rat spots matched to COMMOT: {len(n_r)}/{len(rdf)}")
    n_r.to_parquet(OUT / "rat_spot_table.parquet")

    # combined long-format deliverable
    exo = pd.concat([
        n_h.assign(species="human", treatment=n_h["sample"],
                   timepoint="NA"),
        n_r.assign(species="rat")], ignore_index=True)
    exo.to_parquet(OUT / "exo_scores_per_spot.parquet")
    print(f"[save] exo_scores_per_spot.parquet: {exo.shape}")

    # ============== human statistics (n=1 vs n=1, descriptive) ===========
    hs, ct = n_h[n_h["sample"] == "HS"], n_h[n_h["sample"] == "CTRL"]
    hstats = {
        "framework": ("n=1 vs n=1 proof-of-concept slices; spot-level tests "
                      "descriptive only, no biological replicates"),
        "n_spots": {"HS": int(len(hs)), "CTRL": int(len(ct))},
        "scoring": ("scanpy tl.score_genes ctrl_size=50 on merged object "
                    "(common control-gene pool)"),
        "geneset_provenance": ("7 exosome modules + LDAM/Ferroptosis reference "
                               "sets verbatim from R/27_exosome_function.R "
                               "lines 39-54 (same definitions as analyses "
                               "27/28/36)"),
        "genes_detected": det_h,
        "module_differences": {}, "spearman_vs_gal3": {},
        "spearman_vs_refs": {}, "gal3_receiver_top_quartile": {},
    }
    for m in MODULES + REFS:
        u, p, rbc = rank_biserial(hs[f"score_{m}"].to_numpy(),
                                  ct[f"score_{m}"].to_numpy())
        hstats["module_differences"][m] = {
            "median_HS": float(hs[f"score_{m}"].median()),
            "median_CTRL": float(ct[f"score_{m}"].median()),
            "mannwhitney_U": u, "p_value": p, "rank_biserial_r": rbc}
    # core test: module score x COMMOT Gal-3 sender/receiver per slice
    for tag, sub in [("HS", hs), ("CTRL", ct)]:
        for m in MODULES + REFS:
            for fld in ["gal3_sender_total", "gal3_receiver_total"]:
                r, p, n = spearman(sub[f"score_{m}"].to_numpy(float),
                                   sub[fld].to_numpy(float))
                hstats["spearman_vs_gal3"][f"{m}__{fld}__{tag}"] = {
                    "rho": r, "p_value": p, "n": n}
        for m in MODULES:
            for ref in REFS:
                r, p, n = spearman(sub[f"score_{m}"].to_numpy(float),
                                   sub[f"score_{ref}"].to_numpy(float))
                hstats["spearman_vs_refs"][f"{m}__{ref}__{tag}"] = {
                    "rho": r, "p_value": p, "n": n}
        # co-localisation hotspot: top-quartile Gal-3 receiver spots
        thr = sub["gal3_receiver_total"].quantile(0.75)
        top = sub[sub["gal3_receiver_total"] >= thr]
        rest = sub[sub["gal3_receiver_total"] < thr]
        for m in MODULES + REFS:
            u, p, rbc = rank_biserial(top[f"score_{m}"].to_numpy(),
                                      rest[f"score_{m}"].to_numpy())
            hstats["gal3_receiver_top_quartile"][f"{m}__{tag}"] = {
                "n_top": int(len(top)), "n_rest": int(len(rest)),
                "median_top": float(top[f"score_{m}"].median()),
                "median_rest": float(rest[f"score_{m}"].median()),
                "mannwhitney_U": u, "p_value": p, "rank_biserial_r": rbc}
    json.dump(hstats, open(OUT / "human_stats.json", "w"),
              indent=2, ensure_ascii=False)
    print("[save] human_stats.json")

    # ============== rat time-course + sample-level sync ==================
    tc = json.load(open(ROOT / "results" / "70_commot_spatial" /
                        "timecourse_lr.json"))
    tps = ["D5", "D10", "D20", "D40"]
    rstats = {"framework": ("n=2 per group per timepoint; effect sizes "
                            "primary, spot-level MW descriptive only"),
              "scoring": ("scanpy tl.score_genes ctrl_size=50 per sample; "
                          "human->rat orthologs via ortholog_map.json with "
                          "title-case fallback (provenance recorded)"),
              "per_timepoint": {}, "sample_level": {}}
    for tp in tps:
        se = n_r[(n_r["timepoint"] == tp) & (n_r["treatment"] == "SE")]
        co = n_r[(n_r["timepoint"] == tp) & (n_r["treatment"] == "CTRL")]
        entry = {"n_spots": {"SE": int(len(se)), "CTRL": int(len(co))},
                 "modules": {}}
        for m in MODULES + REFS:
            se_m = se.groupby("sample")[f"score_{m}"].mean()
            co_m = co.groupby("sample")[f"score_{m}"].mean()
            u, p, rbc = rank_biserial(se[f"score_{m}"].to_numpy(),
                                      co[f"score_{m}"].to_numpy())
            entry["modules"][m] = {
                "sample_means_SE": {k: float(v) for k, v in se_m.items()},
                "sample_means_CTRL": {k: float(v) for k, v in co_m.items()},
                "delta_mean_SE_minus_CTRL": float(se_m.mean() - co_m.mean()),
                "spot_level_mannwhitney_p_descriptive": p,
                "spot_level_rank_biserial_r_descriptive": rbc}
        rstats["per_timepoint"][tp] = entry
    json.dump(rstats, open(OUT / "rat_timecourse.json", "w"),
              indent=2, ensure_ascii=False)
    print("[save] rat_timecourse.json")

    # ---- sample-level synchrony (PRIMARY test, n=16) ----
    smean = n_r.groupby("sample").agg(
        treatment=("treatment", "first"), timepoint=("timepoint", "first"),
        **{f"mean_{m}": (f"score_{m}", "mean") for m in MODULES + REFS})
    smean["frac_recv_pos"] = [tc[s]["frac_recv_pos"] for s in smean.index]
    smean["gal3_receiver_mean_pos"] = [tc[s]["gal3_receiver_mean_pos"]
                                       for s in smean.index]
    smean["gal3_receiver_mean"] = [tc[s]["gal3_receiver_mean"]
                                   for s in smean.index]
    sync = {"framework": ("PRIMARY TEST: sample-level Spearman across the 16 "
                          "rat samples (n=16), module mean vs per-sample "
                          "COMMOT Gal-3 broadcast range (frac_recv_pos) and "
                          "intensity (gal3_receiver_mean_pos)"),
            "n_samples": int(len(smean)),
            "per_sample": smean.reset_index().to_dict("records"),
            "spearman_all16": {}, "spearman_SE_only_n8": {}}
    for m in MODULES + REFS:
        for fld in ["frac_recv_pos", "gal3_receiver_mean_pos"]:
            r, p, n = spearman(smean[f"mean_{m}"].to_numpy(float),
                               smean[fld].to_numpy(float))
            sync["spearman_all16"][f"{m}__{fld}"] = {"rho": r, "p": p, "n": n}
            se = smean[smean["treatment"] == "SE"]
            r, p, n = spearman(se[f"mean_{m}"].to_numpy(float),
                               se[fld].to_numpy(float))
            sync["spearman_SE_only_n8"][f"{m}__{fld}"] = {
                "rho": r, "p": p, "n": n}
    # shape comparison: module time-course (SE - CTRL delta) vs COMMOT curve
    shape = {}
    gal3_curve = {}
    for tp in tps:
        e = {}
        for fld in ["frac_recv_pos", "gal3_receiver_mean"]:
            se_v = np.mean([tc[s][fld] for s in tc
                            if tc[s].get("species") == "rat"
                            and tc[s]["treatment"] == "SE"
                            and tc[s]["timepoint"] == tp])
            co_v = np.mean([tc[s][fld] for s in tc
                            if tc[s].get("species") == "rat"
                            and tc[s]["treatment"] == "CTRL"
                            and tc[s]["timepoint"] == tp])
            e[fld] = {"SE": float(se_v), "CTRL": float(co_v),
                      "delta": float(se_v - co_v)}
        gal3_curve[tp] = e
    for m in MODULES:
        delta = {tp: rstats["per_timepoint"][tp]["modules"][m]
                 ["delta_mean_SE_minus_CTRL"] for tp in tps}
        entry = {"module_delta_by_tp": {t: float(delta[t]) for t in tps},
                 "gal3_curves": gal3_curve}
        for fld in ["frac_recv_pos", "gal3_receiver_mean"]:
            r, p = stats.spearmanr(
                [delta[t] for t in tps],
                [gal3_curve[t][fld]["delta"] for t in tps])
            entry[f"spearman_shape_vs_{fld}_n4_descriptive"] = {
                "rho": float(r), "p": float(p)}
        shape[m] = entry
    sync["timecourse_shape_comparison_n4_descriptive"] = shape
    json.dump(sync, open(OUT / "sample_level_sync.json", "w"),
              indent=2, ensure_ascii=False)
    print("[save] sample_level_sync.json")
    print("DONE compute")


# ================================================================== figures
def figures():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    exo = pd.read_parquet(OUT / "exo_scores_per_spot.parquet")
    sync = json.load(open(OUT / "sample_level_sync.json"))
    hs = exo[(exo["species"] == "human") & (exo["sample"] == "HS")]

    def zsc(v):
        return (v - np.nanmean(v)) / np.nanstd(v)

    fig, axes = plt.subplots(4, 2, figsize=(9.0, 15.5))
    rows = [("score_SNARE_fusion", "SNARE_fusion score (z)"),
            ("score_Exosome_uptake", "Exosome_uptake score (z)"),
            ("gal3_receiver_total", "COMMOT Gal-3 receiver signal (z)")]
    panels = [("CTRL", "Human CTRL (GSM8778056)"),
              ("HS", "Human MTLE-HS (GSM8778057)")]
    for i, (col, lab) in enumerate(rows):
        for j, (smp, ttl) in enumerate(panels):
            ax = axes[i, j]
            sub = exo[(exo["species"] == "human") & (exo["sample"] == smp)]
            v = zsc(sub[col].to_numpy(float))
            sc_ = ax.scatter(sub["pxl_col_in_fullres"],
                             -sub["pxl_row_in_fullres"], c=v, s=3,
                             cmap="RdYlBu_r", vmin=-2.5, vmax=2.5,
                             linewidths=0)
            ax.set_title(f"{lab} — {ttl}", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_aspect("equal")
            for sp in ax.spines.values():
                sp.set_visible(False)
        fig.colorbar(sc_, ax=axes[i, :].tolist(), shrink=0.6,
                     label="z-score (shared scale)")
    # bottom: rat sample-level synchrony scatter (n=16)
    smean = pd.DataFrame(sync["per_sample"]).set_index("sample")
    tps = ["D5", "D10", "D20", "D40"]
    mk = {"SE": "o", "CTRL": "s"}
    for j, m in enumerate(["Rab_secretory", "Syntenin_ectosome"]):
        ax = axes[3, j]
        for tp, c in zip(tps, ["#1b9e77", "#d95f02", "#7570b3", "#e7298a"]):
            for tr in ["SE", "CTRL"]:
                sub = smean[(smean["timepoint"] == tp)
                            & (smean["treatment"] == tr)]
                ax.scatter(sub[f"mean_{m}"], sub["frac_recv_pos"], s=45,
                           c=c, marker=mk[tr], edgecolors="k", linewidths=0.4,
                           alpha=0.9,
                           label=f"{tp} {tr}" if j == 0 else None)
        d = sync["spearman_all16"][f"{m}__frac_recv_pos"]
        ax.set_title(f"Rat sample-level (n=16): {m} vs Gal-3 receiver range\n"
                     f"Spearman rho={d['rho']:.2f}, p={d['p']:.3g}",
                     fontsize=9)
        ax.set_xlabel(f"sample mean {m} score")
        ax.set_ylabel("Gal-3 receiver range (frac spots > 0)")
    axes[3, 0].legend(fontsize=6, ncol=2, frameon=False, loc="best")
    fig.suptitle("Exosome machinery x COMMOT Gal-3 broadcast field "
                 "(GSE288888 human; GSE301260 rat)", fontsize=11, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    out = FIGD / "Fig_spatial_exosome.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"[save] {out}")
    print("DONE figs")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--figs-only", action="store_true")
    a = ap.parse_args()
    if not a.figs_only:
        compute()
    figures()
