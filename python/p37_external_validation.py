# -*- coding: utf-8 -*-
"""
p37_external_validation.py
Held-out external validation of the LDAM-ferroptosis axis and the 7-gene
signature in GEO epilepsy cohorts never used in any derivation step.

Cohorts / contrasts:
  GSE186334 (compartment-specific total RNA-seq, postmortem controls):
    1) Hip_cyt  : hippocampus cytoplasm, CTL (n=5) vs mTLE non-HS (n=5)
    2) Cx_cyt_W0: cortex cytoplasm,     CTL (n=6) vs mTLE non-HS (n=6)
    3) Cx_cyt_W4: cortex cytoplasm,     CTL (n=6) vs mTLE+HS (n=6)
  GSE272142 (temporal-lobe RNA-seq):
    4) DG       : TLE dentate gyrus (n=6) vs control temporal lobe (n=6)
  GSE205661 (Agilent-052909 lncRNA+mRNA microarray, GPL19072):
    5) Hip      : TLE+HS hippocampus (n=6) vs normal temporal/parietal
                  cortex from TBI decompression surgery (n=9)

Gene-set provenance:
  - LDAM signature (n=353): genesets/LDAM_geneset.txt
  - Ferroptosis core (n=1301): genesets/Ferroptosis_geneset.txt
  - Module 5 (iron-DAM metacell module, n=28):
    results/59_metacell_wgcna/metacell_modules.csv (module == 5)
  - Hub genes: LGALS3, SPP1, APOE, FTH1
  - 7-gene blood signature: LPCAT2, S100A10, CD9, SYNGR2, ITGAM, ITGB2, YBX1
    (results/04_ML_selected_genes.csv + blood hubs; training cohort GSE143272)

Probe annotation for GPL19072: GEO platform files (GPL19072 SOFT and the
052909 design table) ship with empty gene-annotation columns; the mapping was
recovered from GSE192627_all_expression_signal.xlsx (same Agilent-052909
design, full vendor annotation incl. GeneSymbol), giving 31,860 annotated
probes / 22,621 symbols. Cached in
data/external_validation/GPL19072_probe2symbol.csv.

Statistics: Mann-Whitney U (two-sided), rank-biserial r with bootstrap 95% CI
(5000 resamples, seed 42); logistic model trained on GSE143272 with frozen
coefficients applied to within-cohort z-scored expression; refit apparent AUC
with bootstrap 95% CI.
"""
from pathlib import Path
import json
import gzip
import numpy as np
import pandas as pd
import requests
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究")
DATA = BASE / "data" / "external_validation"
OUT = BASE / "results" / "73_external_validation"
FIGD = BASE / "figures" / "v25"
OUT.mkdir(parents=True, exist_ok=True)
FIGD.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(42)
N_BOOT = 5000

GENES7 = ["LPCAT2", "S100A10", "CD9", "SYNGR2", "ITGAM", "ITGB2", "YBX1"]
HUBS = ["LGALS3", "SPP1", "APOE", "FTH1"]

# ---------------------------------------------------------------- gene sets
def read_geneset(p):
    return [l.strip() for l in open(p, encoding="utf-8") if l.strip()]

LDAM = read_geneset(BASE / "genesets" / "LDAM_geneset.txt")
FERRO = read_geneset(BASE / "genesets" / "Ferroptosis_geneset.txt")
_mm = pd.read_csv(BASE / "results" / "59_metacell_wgcna" / "metacell_modules.csv")
MODULE5 = _mm.loc[_mm["module"] == 5, "gene"].tolist()
GENESETS = {"LDAM": LDAM, "Ferroptosis": FERRO, "Module5": MODULE5}
NEEDED_SYMBOLS = sorted(set(LDAM) | set(FERRO) | set(MODULE5) | set(GENES7) | set(HUBS))
print(f"Gene sets: LDAM={len(LDAM)}, Ferroptosis={len(FERRO)}, Module5={len(MODULE5)}")

# ------------------------------------------------- symbol -> ENSG (mygene)
def map_needed_symbols(symbols, cache):
    if cache.exists():
        d = pd.read_csv(cache)
        return dict(zip(d["ensg"], d["symbol"]))
    sess = requests.Session()
    sess.trust_env = False  # bypass broken system proxy (curl works direct)
    ensg2sym = {}
    syms = list(symbols)
    for i in range(0, len(syms), 60):
        chunk = syms[i:i + 60]
        q = " OR ".join(f'symbol:"{s}"' for s in chunk)
        recs = None
        for attempt in range(3):
            try:
                r = sess.get("https://mygene.info/v3/query",
                             params={"q": q, "fields": "symbol,ensembl.gene",
                                     "species": "human", "size": 1000},
                             timeout=40)
                r.raise_for_status()
                recs = r.json().get("hits", [])
                break
            except Exception as e:
                print(f"  mygene attempt {attempt+1} failed: {type(e).__name__}")
        for hit in recs or []:
            sym = hit.get("symbol")
            ens = hit.get("ensembl")
            ens = ens if isinstance(ens, list) else [ens] if ens else []
            for e in ens:
                if isinstance(e, dict) and e.get("gene"):
                    ensg2sym.setdefault(e["gene"], sym)
        print(f"  symbol-mapped {i + len(chunk)}/{len(syms)}")
    d = pd.DataFrame({"ensg": list(ensg2sym), "symbol": list(ensg2sym.values())})
    d.to_csv(cache, index=False)
    return ensg2sym

# ------------------------------------------------------------- data loaders
def collapse_symbols(df, sym_series):
    """Keep highest-mean row per symbol; returns symbol-indexed df."""
    df = df.copy()
    df["symbol"] = sym_series.values
    df = df.dropna(subset=["symbol"])
    cols = [c for c in df.columns if c != "symbol"]
    df["_mean"] = df[cols].mean(axis=1)
    df = (df.sort_values("_mean", ascending=False)
            .drop_duplicates("symbol")
            .drop(columns="_mean")
            .set_index("symbol"))
    return df

def load_gse186334(fname):
    """Raw counts -> log2 CPM, ENSG->symbol via mygene cache."""
    counts = pd.read_csv(DATA / fname, sep="\t", index_col=0)
    groups = [c.split(".")[0] for c in counts.columns]  # CTL/CTL.1/... -> CTL
    counters, newcols = {}, []
    for g in groups:
        counters[g] = counters.get(g, 0) + 1
        newcols.append(f"{g}_{counters[g]}")
    counts.columns = newcols
    mapping = map_needed_symbols(NEEDED_SYMBOLS, OUT / "ensg2symbol.csv")
    lib = counts.sum(axis=0)
    logcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)
    expr = collapse_symbols(logcpm, logcpm.index.to_series().map(mapping))
    return expr, np.array(groups)

def load_gse272142():
    """Annotated raw counts xlsx (converted to csv); Symbol column."""
    df = pd.read_csv(DATA / "GSE272142_genes.rawcount.anno.csv")
    samples = ["TLE1_DG", "TLE3_DG", "TLE5_DG", "TLE6_DG", "TLE7_DG", "TLE8_DG",
               "XY161222", "XY171205", "XY171229", "XY180413", "XYA0513A",
               "XYA0516B"]
    counts = df.set_index("#id")[samples]
    lib = counts.sum(axis=0)
    logcpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)
    expr = collapse_symbols(logcpm, df["Symbol"])
    groups = np.array(["TLE"] * 6 + ["CTL"] * 6)
    return expr, groups

def load_gse205661():
    """GPL19072 series matrix (GeneSpring quantile-normalised log2 values),
    probe->symbol via recovered vendor annotation (GSE192627)."""
    p2s = pd.read_csv(DATA / "GPL19072_probe2symbol.csv")
    p2s = p2s.dropna(subset=["symbol"])
    probe2sym = dict(zip(p2s["probe"], p2s["symbol"]))
    rows, header = [], None
    with gzip.open(DATA / "GSE205661-GPL19072_series_matrix.txt.gz", "rt",
                   encoding="utf-8", errors="replace") as fh:
        in_table = False
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            parts = line.rstrip("\n").split("\t")
            parts = [p.strip('"') for p in parts]
            if header is None:
                header = parts
                continue
            rows.append(parts)
    mat = pd.DataFrame(rows, columns=header).set_index("ID_REF")
    mat = mat.apply(pd.to_numeric, errors="coerce")
    # sample titles: first 6 = human_RNA_TLE+HS, next 9 = human_RNA_NTP
    groups = np.array(["TLE" if i < 6 else "CTL" for i in range(mat.shape[1])])
    expr = collapse_symbols(mat, mat.index.to_series().map(probe2sym))
    return expr, groups

CONTRASTS = {
    "GSE186334_Hip_cyt_CTL_vs_mTLEnonHS": {
        "loader": lambda: load_gse186334("GSE186334_Hip_cyt_CTLvsW0_countTable_edgeR.txt.gz"),
        "case_label": "W0", "tissue": "Hippocampus (cytoplasm)",
        "case_name": "mTLE non-HS", "ctrl_name": "postmortem control",
        "platform": "Illumina total RNA-seq"},
    "GSE186334_Cx_cyt_CTL_vs_mTLEnonHS": {
        "loader": lambda: load_gse186334("GSE186334_Cx_cyt_CTLvsW0_countTable_edgeR.txt.gz"),
        "case_label": "W0", "tissue": "Cortex (cytoplasm)",
        "case_name": "mTLE non-HS", "ctrl_name": "postmortem control",
        "platform": "Illumina total RNA-seq"},
    "GSE186334_Cx_cyt_CTL_vs_mTLEHS": {
        "loader": lambda: load_gse186334("GSE186334_Cx_cyt_CTLvsW4_countTable_edgeR.txt.gz"),
        "case_label": "W4", "tissue": "Cortex (cytoplasm)",
        "case_name": "mTLE+HS", "ctrl_name": "postmortem control",
        "platform": "Illumina total RNA-seq"},
    "GSE272142_DG_CTL_vs_TLE": {
        "loader": load_gse272142,
        "case_label": "TLE", "tissue": "Dentate gyrus (temporal lobe)",
        "case_name": "TLE", "ctrl_name": "control temporal lobe",
        "platform": "Illumina RNA-seq"},
    "GSE205661_Hip_NTP_vs_TLEHS": {
        "loader": load_gse205661,
        "case_label": "TLE", "tissue": "Hippocampus vs temporal/parietal cortex",
        "case_name": "TLE+HS", "ctrl_name": "NTP (TBI decompression)",
        "platform": "Agilent-052909 microarray (GPL19072)"},
}

# ------------------------------------------------------------------ stats
def rank_biserial(x_case, x_ctrl):
    n1, n0 = len(x_case), len(x_ctrl)
    U = mannwhitneyu(x_case, x_ctrl, alternative="two-sided").statistic
    return 2 * U / (n1 * n0) - 1

def boot_ci(case, ctrl, fun, n=N_BOOT, rng=RNG):
    stats = []
    case, ctrl = np.asarray(case), np.asarray(ctrl)
    for _ in range(n):
        b1 = rng.choice(case, size=len(case), replace=True)
        b0 = rng.choice(ctrl, size=len(ctrl), replace=True)
        try:
            stats.append(fun(b1, b0))
        except Exception:
            pass
    lo, hi = np.percentile(stats, [2.5, 97.5])
    return float(lo), float(hi)

def mw_summary(case, ctrl):
    p = float(mannwhitneyu(case, ctrl, alternative="two-sided").pvalue)
    r = float(rank_biserial(case, ctrl))
    lo, hi = boot_ci(case, ctrl, lambda a, b: rank_biserial(a, b))
    return {"p_mw": p, "rank_biserial_r": r, "r_ci95": [lo, hi],
            "mean_case": float(np.mean(case)), "mean_ctrl": float(np.mean(ctrl))}

def auc_ci(y, score):
    a = float(roc_auc_score(y, score))
    idx = np.arange(len(y))
    stats = []
    for _ in range(N_BOOT):
        b = RNG.choice(idx, size=len(idx), replace=True)
        if len(set(y[b])) < 2:
            continue
        stats.append(roc_auc_score(y[b], score[b]))
    lo, hi = np.percentile(stats, [2.5, 97.5])
    return a, float(lo), float(hi)

# -------------------------------------------------- train frozen 7-gene model
train = pd.read_csv(BASE / "results" / "08_GSE143272_normalized.csv", index_col=0)
meta = pd.read_csv(BASE / "results" / "08_GSE143272_sample_meta.csv")
meta = meta.set_index("sample").loc[train.columns]
y_train = (meta["group"] == "Epilepsy").astype(int).values
Xtr = train.loc[GENES7].T.values.astype(float)
mu, sd = Xtr.mean(axis=0), Xtr.std(axis=0, ddof=1)
Ztr = (Xtr - mu) / sd
clf = LogisticRegression(C=1e6, max_iter=5000)
clf.fit(Ztr, y_train)
coefs = dict(zip(GENES7, clf.coef_[0].tolist()))
intercept = float(clf.intercept_[0])
train_auc = float(roc_auc_score(y_train, clf.decision_function(Ztr)))
print(f"Frozen model trained on GSE143272: apparent AUC={train_auc:.3f}")
print("coefficients:", {k: round(v, 3) for k, v in coefs.items()})

# ------------------------------------------------------------------ main loop
all_results = {"training": {"cohort": "GSE143272",
                            "n_epilepsy": int(y_train.sum()),
                            "n_control": int((1 - y_train).sum()),
                            "frozen_coefs": coefs, "intercept": intercept,
                            "apparent_auc_train": train_auc},
               "gene_sets": {"LDAM": {"source": "genesets/LDAM_geneset.txt", "n": len(LDAM)},
                             "Ferroptosis": {"source": "genesets/Ferroptosis_geneset.txt", "n": len(FERRO)},
                             "Module5": {"source": "results/59_metacell_wgcna/metacell_modules.csv (module==5)",
                                         "n": len(MODULE5)},
                             "signature_7gene": GENES7, "hubs": HUBS},
               "contrasts": {}}

roc_data_frozen, roc_data_refit = {}, {}

for cname, cfg in CONTRASTS.items():
    print(f"\n=== {cname} ===")
    expr, groups = cfg["loader"]()
    print(f"  matrix: {expr.shape[0]} symbols x {expr.shape[1]} samples")
    case_idx = np.where(groups == cfg["case_label"])[0]
    ctrl_idx = np.where(groups == "CTL")[0]
    assert len(case_idx) > 0 and len(ctrl_idx) > 0, "group parsing failed"
    y = np.zeros(expr.shape[1], dtype=int)
    y[case_idx] = 1

    res = {"design": {"tissue": cfg["tissue"], "case": cfg["case_name"],
                      "ctrl": cfg["ctrl_name"], "platform": cfg["platform"],
                      "n_case": int(len(case_idx)), "n_ctrl": int(len(ctrl_idx))},
           "axis_scores": {}, "hubs": {}, "signature": {}}

    Z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1, ddof=1), axis=0)
    Z = Z.dropna()  # drop constant genes (std = 0)

    # (a) axis scores
    for sname, genes in GENESETS.items():
        det = [g for g in genes if g in Z.index]
        score = Z.loc[det].mean(axis=0).values
        s_case, s_ctrl = score[case_idx], score[ctrl_idx]
        summ = mw_summary(s_case, s_ctrl)
        summ["n_detected"] = len(det)
        summ["coverage"] = float(len(det) / len(genes))
        summ["scores"] = {"case": s_case.tolist(), "ctrl": s_ctrl.tolist()}
        res["axis_scores"][sname] = summ
        print(f"  {sname}: r={summ['rank_biserial_r']:+.3f} "
              f"[{summ['r_ci95'][0]:+.3f},{summ['r_ci95'][1]:+.3f}] "
              f"P={summ['p_mw']:.4f} (genes {len(det)}/{len(genes)})")

    # (c) hub genes
    for g in HUBS:
        if g in expr.index:
            v = expr.loc[g].values
            summ = mw_summary(v[case_idx], v[ctrl_idx])
            summ["log2FC"] = float(v[case_idx].mean() - v[ctrl_idx].mean())
            res["hubs"][g] = summ
            print(f"  hub {g}: log2FC={summ['log2FC']:+.3f} P={summ['p_mw']:.4f}")
        else:
            res["hubs"][g] = None
            print(f"  hub {g}: NOT DETECTED")

    # (b) 7-gene signature
    det7 = [g for g in GENES7 if g in Z.index]
    res["signature"]["genes_detected"] = det7
    if len(det7) == len(GENES7):
        Zext = Z.loc[GENES7].T.values
        frozen_score = Zext @ clf.coef_[0]
        a_f, lo_f, hi_f = auc_ci(y, frozen_score)
        res["signature"]["frozen_coef_model"] = {
            "auc": a_f, "auc_ci95": [lo_f, hi_f],
            "note": "GSE143272-trained coefficients applied to within-cohort "
                    "z-scored expression; intercept omitted (ranking only)."}
        clf2 = LogisticRegression(C=1e6, max_iter=5000)
        clf2.fit(Zext, y)
        refit_score = clf2.decision_function(Zext)
        a_r, lo_r, hi_r = auc_ci(y, refit_score)
        res["signature"]["refit_model"] = {
            "auc": a_r, "auc_ci95": [lo_r, hi_r],
            "coefs": dict(zip(GENES7, clf2.coef_[0].tolist()))}
        print(f"  7-gene: frozen AUC={a_f:.3f} [{lo_f:.3f},{hi_f:.3f}]  "
              f"refit AUC={a_r:.3f} [{lo_r:.3f},{hi_r:.3f}]")
        fpr, tpr, _ = roc_curve(y, frozen_score)
        roc_data_frozen[cname] = (fpr, tpr, a_f, lo_f, hi_f)
        fpr2, tpr2, _ = roc_curve(y, refit_score)
        roc_data_refit[cname] = (fpr2, tpr2, a_r, lo_r, hi_r)
        sig_dir = {}
        for g in GENES7:
            v = expr.loc[g].values
            sig_dir[g] = {"log2FC": float(v[case_idx].mean() - v[ctrl_idx].mean()),
                          "p_mw": float(mannwhitneyu(v[case_idx], v[ctrl_idx],
                                                     alternative="two-sided").pvalue)}
        res["signature"]["gene_level"] = sig_dir
    else:
        res["signature"]["error"] = f"only {len(det7)}/7 genes detected"
        print(f"  7-gene: MISSING genes -> {set(GENES7) - set(det7)}")

    all_results["contrasts"][cname] = res
    sc = pd.DataFrame({"sample": expr.columns,
                       "group": np.where(y == 1, "TLE", "CTL")})
    for sname in GENESETS:
        det = [g for g in GENESETS[sname] if g in Z.index]
        sc[sname + "_score"] = Z.loc[det].mean(axis=0).values
    sc.to_csv(OUT / f"73_{cname}_axis_scores.csv", index=False)

# -------------------------------------------------- direction consistency
prior_axis = {"LDAM": +1, "Ferroptosis": +1, "Module5": +1}
prior_hub = {"LGALS3": +1, "SPP1": +1, "APOE": +1, "FTH1": +1}     # GSE190451
prior_hub_gse168375 = {"LGALS3": +1, "SPP1": +1, "APOE": +1, "FTH1": -1}
consistency = {"axis": [], "hubs_vs_GSE190451": [], "hubs_vs_GSE168375": [],
               "signature_auc_gt_0.5": []}
for cname, res in all_results["contrasts"].items():
    for sname in GENESETS:
        r = res["axis_scores"][sname]["rank_biserial_r"]
        consistency["axis"].append(
            {"contrast": cname, "score": sname, "r": r,
             "prior_direction": prior_axis[sname],
             "concordant": bool(np.sign(r) == prior_axis[sname])})
    for g in HUBS:
        h = res["hubs"][g]
        if h:
            for key, prior in (("hubs_vs_GSE190451", prior_hub),
                               ("hubs_vs_GSE168375", prior_hub_gse168375)):
                consistency[key].append(
                    {"contrast": cname, "gene": g, "log2FC": h["log2FC"],
                     "concordant": bool(np.sign(h["log2FC"]) == prior[g])})
    sig = res["signature"].get("frozen_coef_model")
    if sig:
        consistency["signature_auc_gt_0.5"].append(
            {"contrast": cname, "auc_frozen": sig["auc"],
             "ci_excludes_0.5": sig["auc_ci95"][0] > 0.5,
             "concordant": sig["auc"] > 0.5})
all_results["direction_consistency"] = consistency

# -------------------------------------------------- cross-cohort summary table
rows = []
for cname, res in all_results["contrasts"].items():
    row = {"contrast": cname,
           "n_case": res["design"]["n_case"], "n_ctrl": res["design"]["n_ctrl"]}
    for sname in GENESETS:
        s = res["axis_scores"][sname]
        row[f"{sname}_r"] = round(s["rank_biserial_r"], 3)
        row[f"{sname}_P"] = round(s["p_mw"], 4)
    sig = res["signature"].get("frozen_coef_model", {})
    row["AUC_frozen"] = round(sig.get("auc", np.nan), 3)
    row["AUC_CI"] = (f"{sig['auc_ci95'][0]:.2f}-{sig['auc_ci95'][1]:.2f}"
                     if sig else "NA")
    hc = [c for c in consistency["hubs_vs_GSE190451"] if c["contrast"] == cname]
    row["hub_dir_concordance"] = f"{sum(c['concordant'] for c in hc)}/{len(hc)}"
    rows.append(row)
summary_df = pd.DataFrame(rows)
summary_df.to_csv(OUT / "73_cross_cohort_summary.csv", index=False)
all_results["cross_cohort_summary"] = rows

n_conc = sum(c["concordant"] for c in consistency["axis"])
all_results["direction_consistency_summary"] = {
    "axis_concordance": f"{n_conc}/{len(consistency['axis'])}",
    "hub_concordance_vs_GSE190451":
        f"{sum(c['concordant'] for c in consistency['hubs_vs_GSE190451'])}/{len(consistency['hubs_vs_GSE190451'])}",
    "hub_concordance_vs_GSE168375":
        f"{sum(c['concordant'] for c in consistency['hubs_vs_GSE168375'])}/{len(consistency['hubs_vs_GSE168375'])}",
    "signature_auc_above_chance":
        f"{sum(c['concordant'] for c in consistency['signature_auc_gt_0.5'])}/{len(consistency['signature_auc_gt_0.5'])}",
    "signature_auc_ci_excludes_chance":
        f"{sum(c['ci_excludes_0.5'] for c in consistency['signature_auc_gt_0.5'])}/{len(consistency['signature_auc_gt_0.5'])}",
}
print("\nDirection consistency:",
      json.dumps(all_results["direction_consistency_summary"], indent=2))

with open(OUT / "validation_stats.json", "w", encoding="utf-8") as f:
    json.dump(all_results, f, indent=2, ensure_ascii=False)
print(f"Stats saved: {OUT / 'validation_stats.json'}")

# ------------------------------------------------------------------ figure
SHORT = {
    "GSE186334_Hip_cyt_CTL_vs_mTLEnonHS": ("GSE186334", "Hippocampus: CTL vs mTLE non-HS"),
    "GSE186334_Cx_cyt_CTL_vs_mTLEnonHS": ("GSE186334", "Cortex: CTL vs mTLE non-HS"),
    "GSE186334_Cx_cyt_CTL_vs_mTLEHS": ("GSE186334", "Cortex: CTL vs mTLE+HS"),
    "GSE272142_DG_CTL_vs_TLE": ("GSE272142", "Dentate gyrus: CTL vs TLE"),
    "GSE205661_Hip_NTP_vs_TLEHS": ("GSE205661", "Hippocampus: NTP vs TLE+HS"),
}
NROWS = len(CONTRASTS)
fig = plt.figure(figsize=(13, 4 + 3.1 * NROWS))
gs = fig.add_gridspec(2, 1, height_ratios=[NROWS, 2.1], hspace=0.30)
sub1 = fig.add_subfigure(gs[0])
axs = sub1.subplots(NROWS, 3)
sub1.subplots_adjust(hspace=0.95, wspace=0.30)
for row, (cname, cfg) in enumerate(CONTRASTS.items()):
    res = all_results["contrasts"][cname]
    for si, sname in enumerate(GENESETS):
        summ = res["axis_scores"][sname]
        s_case = np.array(summ["scores"]["case"])
        s_ctrl = np.array(summ["scores"]["ctrl"])
        ax = axs[row, si]
        bp = ax.boxplot([s_ctrl, s_case], widths=0.55, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black", lw=1.4))
        for patch, c in zip(bp["boxes"], ["#9ecae1", "#fb9a99"]):
            patch.set_facecolor(c)
        for j, (vals, c) in enumerate(zip([s_ctrl, s_case], ["#08519c", "#cb181d"])):
            ax.scatter(np.full(len(vals), j + 1) + RNG.normal(0, 0.05, len(vals)),
                       vals, color=c, s=24, zorder=3)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"CTL (n={len(s_ctrl)})", f"TLE (n={len(s_case)})"],
                           fontsize=8)
        dset, short = SHORT[cname]
        ax.set_title(f"{sname} | {dset} {short}\n"
                     f"r = {summ['rank_biserial_r']:+.2f} "
                     f"[{summ['r_ci95'][0]:+.2f}, {summ['r_ci95'][1]:+.2f}], "
                     f"P = {summ['p_mw']:.3g}", fontsize=8.3)
        if si == 0:
            ax.set_ylabel("Axis score (mean z)", fontsize=8.5)
sub1.suptitle("A  Axis score distributions in held-out cohorts (GSE186334, "
              "GSE272142, GSE205661)", fontsize=12, fontweight="bold",
              x=0.02, ha="left")
sub2 = fig.add_subfigure(gs[1])
axr2 = sub2.subplots(1, 2)
sub2.subplots_adjust(wspace=0.25)
for ax, rd, ttl in [(axr2[0], roc_data_frozen, "frozen GSE143272 coefficients"),
                    (axr2[1], roc_data_refit, "refit within external cohort")]:
    for cname, (fpr, tpr, a, lo, hi) in rd.items():
        dset, short = SHORT[cname]
        ax.plot(fpr, tpr, lw=1.6,
                label=f"{dset} {short}: AUC = {a:.2f} ({lo:.2f}-{hi:.2f})")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(f"7-gene signature ROC — {ttl}", fontsize=10)
    ax.legend(fontsize=6.8, loc="lower right")
sub2.suptitle("B  Seven-gene signature separability (LPCAT2, S100A10, CD9, "
              "SYNGR2, ITGAM, ITGB2, YBX1)", fontsize=12, fontweight="bold",
              x=0.02, ha="left")
figpath = FIGD / "Fig_external_validation.png"
fig.savefig(figpath, dpi=300, bbox_inches="tight")
plt.close("all")
print(f"Figure saved: {figpath}")
