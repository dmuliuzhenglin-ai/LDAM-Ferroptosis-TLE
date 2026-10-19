# -*- coding: utf-8 -*-
"""
p38_sex_stratified.py
Sex-stratified sensitivity analysis (sex as a biological variable) across bulk cohorts.

Cohorts:
  GSE134697  TLE hippocampus vs neocortex (17 paired patients)   [sex inferred: XIST/UTY]
  GSE168375  epilepsy vs control (hippocampus RNA-seq)           [sex inferred: XIST/UTY]
  GSE190451  3 TLE vs 3 control (composition only)               [sex inferred: XIST/UTY]
  GSE63808   129 TLE hippocampus microarray, no controls         [sex inferred: XIST/UTY probes]
  GSE143272  blood, 34 epilepsy vs 51 control                    [sex from GEO metadata]

Gene-set provenance:
  LDAM        : genesets/LDAM_geneset.txt (n=353)
  Ferroptosis : genesets/Ferroptosis_geneset.txt (n=1301, FerrDb-derived)
  WGCNA brown : results/03_WGCNA_LDAM_hub_module_genes.txt (n=582, LDAM hub module MEbrown,
                see results/03_WGCNA_module_trait_correlation.csv)
  7-gene signature: LPCAT2,S100A10,CD9,SYNGR2,ITGAM,ITGB2,YBX1
                (results/04_ML_selected_genes.csv + results/33_blood_exosome/33_blood_hub_genes.csv,
                 combined logit in results/37_blood_combined_ROC; scripts/p34_xgb_shap.py)
  Hub genes   : LGALS3, SPP1, APOE, FTH1

Outputs: results/74_sex_stratified/{sex_composition.csv, stratified_stats.json,
         interaction_tests.json, report.md}, figures/v25/Fig_sex_stratified.png
"""
import json
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import statsmodels.formula.api as smf

BASE = Path(r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究")
OUT = BASE / "results" / "74_sex_stratified"
FIG = BASE / "figures" / "v25"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

MIN_TEST_N = 6   # per-class minimum for MW/AUC
MIN_INTERACTION_N = 8  # per-sex stratum minimum for interaction model

Y_GENES = ["UTY", "RPS4Y1", "KDM5D", "DDX3Y", "EIF1AY", "USP9Y", "ZFY", "NLGN4Y"]
HUB_GENES = ["LGALS3", "SPP1", "APOE", "FTH1"]
GENES7 = ["LPCAT2", "S100A10", "CD9", "SYNGR2", "ITGAM", "ITGB2", "YBX1"]

LDAM_SET = pd.read_csv(BASE / "genesets" / "LDAM_geneset.txt", header=None)[0].str.strip().tolist()
FERRO_SET = pd.read_csv(BASE / "genesets" / "Ferroptosis_geneset.txt", header=None)[0].str.strip().tolist()
MODULE_SET = pd.read_csv(BASE / "results" / "03_WGCNA_LDAM_hub_module_genes.txt", header=None)[0].str.strip().tolist()


def load_expr(path):
    e = pd.read_csv(path)
    e = e.rename(columns={e.columns[0]: "gene"})
    e["gene"] = e["gene"].astype(str).str.strip()
    return e.groupby("gene").mean(numeric_only=True)


def zscore_rows(mat):
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0)
    return z.fillna(0.0)


def geneset_score(expr_z, genes):
    g = [x for x in genes if x in expr_z.index]
    if len(g) < 5:
        return None, len(g)
    return expr_z.loc[g].mean(axis=0), len(g)


def infer_sex_rnaseq(expr):
    """Infer sex from XIST vs Y-linked gene expression. Returns (Series F/M, margin Series)."""
    avail_y = [g for g in Y_GENES if g in expr.index]
    if "XIST" not in expr.index or len(avail_y) == 0:
        return None, None, avail_y
    z = zscore_rows(expr.loc[["XIST"] + avail_y])
    xist = z.loc["XIST"]
    ymean = z.loc[avail_y].mean(axis=0)
    margin = xist - ymean
    sex = pd.Series(np.where(margin > 0, "Female", "Male"), index=expr.columns)
    return sex, margin, avail_y


def rank_biserial_and_mw(case, ctrl):
    """Mann-Whitney + rank-biserial r (positive = higher in case)."""
    U, p = stats.mannwhitneyu(case, ctrl, alternative="two-sided")
    r = 2.0 * U / (len(case) * len(ctrl)) - 1.0  # U computed with case first
    return float(U), float(p), float(r)


# ---------------- load cohorts ----------------
cohorts = {}  # name -> dict(df with group/sex/scores..., sex_source, notes)

# GSE134697
e134 = load_expr(BASE / "results" / "01_GSE134697_vst_normalized.csv")
m134 = pd.read_csv(BASE / "results" / "01_GSE134697_sample_meta.csv")
e134.columns = [c.replace(".", "-") for c in e134.columns]
m134["group"] = np.where(m134["tissue"] == "Hippocampus", "TLE", "Control")
sex134, margin134, y134 = infer_sex_rnaseq(e134)
cohorts["GSE134697"] = dict(
    expr=e134, meta=m134.set_index("sample"), sex=sex134, margin=margin134,
    sex_source="inferred (XIST/UTY)",
    note="paired hippocampus(TLE) vs neocortex(control); MW treats pairs as independent (approximation)",
    case_label="TLE", ctrl_label="Control",
)

# GSE168375
e168 = load_expr(BASE / "results" / "02_GSE168375_vst_normalized.csv")
m168 = pd.read_csv(BASE / "results" / "02_GSE168375_sample_meta.csv")
m168["group"] = np.where(m168["group"].str.lower() == "epilepsy", "TLE", "Control")
sex168, margin168, y168 = infer_sex_rnaseq(e168)
cohorts["GSE168375"] = dict(
    expr=e168, meta=m168.set_index("sample"), sex=sex168, margin=margin168,
    sex_source="inferred (XIST/UTY)", note="hippocampus bulk RNA-seq",
    case_label="TLE", ctrl_label="Control",
)

# GSE190451 (3v3 - composition only)
e190 = load_expr(BASE / "results" / "02_GSE190451_normalized.csv")
m190 = pd.read_csv(BASE / "results" / "02_GSE190451_sample_meta.csv")
m190["group"] = np.where(m190["group"].str.contains("TLE"), "TLE", "Control")
sex190, margin190, y190 = infer_sex_rnaseq(e190)
cohorts["GSE190451"] = dict(
    expr=e190, meta=m190.set_index("sample"), sex=sex190, margin=margin190,
    sex_source="inferred (XIST/UTY)", note="only 3 TLE vs 3 control: composition reported, no tests",
    case_label="TLE", ctrl_label="Control",
)

# GSE143272 blood - real sex
e143 = load_expr(BASE / "results" / "08_GSE143272_normalized.csv")
m143 = pd.read_csv(BASE / "results" / "08_GSE143272_sample_meta.csv")
m143["group"] = np.where(m143["group"].str.lower() == "epilepsy", "TLE", "Control")
g143sex = pd.read_csv(BASE / "data" / "sex_meta" / "GSE143272_sex.csv")


def norm_title(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower().encode("ascii", "ignore").decode())


g143sex["key"] = g143sex["title"].map(norm_title)
sex_map = g143sex.set_index("key")["sex"].to_dict()
m143["sex"] = [sex_map.get(norm_title(s), None) for s in m143["sample"]]
cohorts["GSE143272"] = dict(
    expr=e143, meta=m143.set_index("sample"), sex=None,
    sex_source="GEO metadata", note="whole blood microarray",
    case_label="TLE", ctrl_label="Control",
)

# GSE63808 - probe-level, all TLE, sex inferred, composition + sex-effect only
gsm638 = pd.read_csv(BASE / "data" / "sex_meta" / "GSE63808_pheno.csv")
scores638 = pd.read_csv(BASE / "results" / "03_GSE63808_ssGSEA_scores.csv")
annot = {}
with gzip.open(BASE / "GSE63808" / "GPL6947.annot.gz", "rt", errors="ignore") as f:
    header = None
    for ln in f:
        if ln.startswith("ID\t"):
            header = ln.rstrip("\n").split("\t")
            continue
        if header and not ln.startswith("!"):
            parts = ln.rstrip("\n").split("\t")
            annot[parts[0]] = parts[2] if len(parts) > 2 else ""
        if ln.startswith("!platform_table_end"):
            break
xist_probes = [p for p, s in annot.items() if s == "XIST"]
y_probes = [p for p, s in annot.items() if s in Y_GENES]
sex_rows = {}
with gzip.open(BASE / "GSE63808" / "GSE63808_series_matrix.txt.gz", "rt", errors="ignore") as f:
    table = False
    cols = None
    for ln in f:
        if ln.startswith("!series_matrix_table_begin"):
            table = True
            continue
        if ln.startswith("!series_matrix_table_end"):
            break
        if table:
            parts = ln.rstrip("\n").split("\t")
            if cols is None:
                cols = [c.strip('"') for c in parts[1:]]
                continue
            pid = parts[0].strip('"')
            if pid in xist_probes or pid in y_probes:
                vals = np.array([float(x.strip('"')) for x in parts[1:]])
                sex_rows[annot[pid]] = vals  # last probe per symbol wins
sx = pd.DataFrame(sex_rows, index=cols)
sx = np.log2(sx + 1.0)
zx = (sx - sx.mean()) / sx.std()
xist638 = zx["XIST"] if "XIST" in zx.columns else None
ygenes638 = [g for g in Y_GENES if g in zx.columns]
ymean638 = zx[ygenes638].mean(axis=1)
margin638 = xist638 - ymean638
sex638 = pd.Series(np.where(margin638 > 0, "Female", "Male"), index=cols)
meta638 = pd.DataFrame({"sample": cols, "group": "TLE", "sex": sex638.values}).set_index("sample")
meta638 = meta638.join(scores638.set_index("sample")[["LDAM", "Ferroptosis"]], how="left")
cohorts["GSE63808"] = dict(
    expr=None, meta=meta638, sex=None,
    sex_source="inferred (XIST/UTY probes)",
    note="129 TLE hippocampus, NO control group -> composition + between-sex descriptive only",
    case_label="TLE", ctrl_label=None,
)

# ---------------- composition tables ----------------
comp_rows = []
for name, c in cohorts.items():
    meta = c["meta"].copy()
    if "sex" not in meta.columns or meta["sex"].isna().all():
        meta["sex"] = c["sex"].reindex(meta.index).values if c["sex"] is not None else None
    tab = meta.dropna(subset=["sex"]).groupby(["group", "sex"]).size().unstack(fill_value=0)
    for col in ["Female", "Male"]:
        if col not in tab.columns:
            tab[col] = 0
    fisher_p = np.nan
    if c["ctrl_label"] is not None and tab.shape[0] == 2 and tab.sum().sum() > 0:
        try:
            _, fisher_p = stats.fisher_exact(
                tab.loc[[c["case_label"], c["ctrl_label"]], ["Female", "Male"]].values)
        except Exception:
            pass
    comp_rows.append(dict(
        cohort=name,
        sex_source=c["sex_source"],
        case_F=int(tab.loc[c["case_label"], "Female"]) if c["case_label"] in tab.index else 0,
        case_M=int(tab.loc[c["case_label"], "Male"]) if c["case_label"] in tab.index else 0,
        ctrl_F=int(tab.loc[c["ctrl_label"], "Female"]) if c["ctrl_label"] and c["ctrl_label"] in tab.index else np.nan,
        ctrl_M=int(tab.loc[c["ctrl_label"], "Male"]) if c["ctrl_label"] and c["ctrl_label"] in tab.index else np.nan,
        total=int(tab.sum().sum()),
        fisher_sex_imbalance_p=fisher_p,
        note=c["note"],
    ))
comp = pd.DataFrame(comp_rows)
comp.to_csv(OUT / "sex_composition.csv", index=False)
print(comp.to_string(index=False))

# QC of sex inference: margin summaries
for name in ["GSE134697", "GSE168375", "GSE190451"]:
    m = cohorts[name]["margin"]
    if m is not None:
        f_m = m[cohorts[name]["sex"] == "Female"]
        m_m = m[cohorts[name]["sex"] == "Male"]
        print(f"{name} inference QC: Female margin {f_m.min():.2f}..{f_m.max():.2f}; "
              f"Male margin {m_m.min():.2f}..{m_m.max():.2f}")
f638 = margin638[sex638 == "Female"]; m638 = margin638[sex638 == "Male"]
print(f"GSE63808 inference QC: Female margin {f638.min():.2f}..{f638.max():.2f}; "
      f"Male margin {m638.min():.2f}..{m638.max():.2f}")

# ---------------- per-cohort feature table ----------------
strat_stats = {}
interaction_tests = {}
plot_data = []

for name, c in cohorts.items():
    meta = c["meta"].copy()
    if "sex" not in meta.columns or meta["sex"].isna().all():
        meta["sex"] = c["sex"].reindex(meta.index).values if c["sex"] is not None else None
    meta = meta.dropna(subset=["sex"])

    feats = {}
    if c["expr"] is not None:
        ez = zscore_rows(c["expr"])
        for setname, genes in [("LDAM", LDAM_SET), ("Ferroptosis", FERRO_SET),
                               ("WGCNA_brown_module", MODULE_SET)]:
            sc, ng = geneset_score(ez, genes)
            if sc is not None:
                feats[setname] = sc
                print(f"{name} {setname}: {ng} genes used")
        for g in HUB_GENES:
            if g in c["expr"].index:
                feats[g] = c["expr"].loc[g]
        # 7-gene signature: logistic fit on full cohort (z-scored genes), per-sample score
        avail7 = [g for g in GENES7 if g in ez.index]
        if len(avail7) == 7 and c["ctrl_label"] is not None:
            X7 = ez.loc[GENES7].T
            y7 = (meta.loc[X7.index, "group"] == c["case_label"]).astype(int).values
            if len(np.unique(y7)) == 2:
                lr = LogisticRegression(max_iter=1000, C=1.0)
                lr.fit(X7.values, y7)
                feats["Signature7_logit"] = pd.Series(lr.decision_function(X7.values), index=X7.index)
    else:
        # GSE63808: precomputed ssGSEA scores only
        for s in ["LDAM", "Ferroptosis"]:
            if s in meta.columns:
                feats[s] = meta[s]

    fdf = pd.DataFrame(feats)
    fdf = fdf.join(meta[["group", "sex"]], how="inner")

    # ---- stratified stats ----
    cres = {}
    if c["ctrl_label"] is not None:
        for sexv in ["Female", "Male"]:
            sub = fdf[fdf["sex"] == sexv]
            case = sub[sub["group"] == c["case_label"]]
            ctrl = sub[sub["group"] == c["ctrl_label"]]
            sres = {"n_case": int(len(case)), "n_ctrl": int(len(ctrl))}
            for feat in fdf.columns:
                if feat in ["group", "sex"]:
                    continue
                ca, ct = case[feat].dropna(), ctrl[feat].dropna()
                entry = {"n_case": int(len(ca)), "n_ctrl": int(len(ct)),
                         "median_case": float(ca.median()) if len(ca) else None,
                         "median_ctrl": float(ct.median()) if len(ct) else None}
                if len(ca) >= MIN_TEST_N and len(ct) >= MIN_TEST_N:
                    U, p, r = rank_biserial_and_mw(ca.values, ct.values)
                    entry.update(mw_U=U, mw_p=p, rank_biserial_r=r)
                    entry["AUC"] = float(roc_auc_score(
                        [1] * len(ca) + [0] * len(ct), list(ca.values) + list(ct.values)))
                else:
                    entry["test"] = "skipped (n<6 per class): descriptive only"
                sres[feat] = entry
            cres[sexv] = sres
    else:
        # no control: between-sex comparison of scores (descriptive + MW)
        for feat in fdf.columns:
            if feat in ["group", "sex"]:
                continue
            fv = fdf.loc[fdf["sex"] == "Female", feat].dropna()
            mv = fdf.loc[fdf["sex"] == "Male", feat].dropna()
            entry = {"n_F": int(len(fv)), "n_M": int(len(mv)),
                     "median_F": float(fv.median()), "median_M": float(mv.median())}
            if len(fv) >= MIN_TEST_N and len(mv) >= MIN_TEST_N:
                U, p, r = rank_biserial_and_mw(mv.values, fv.values)
                entry.update(mw_U=U, mw_p=p, rank_biserial_r_MvsF=r)
            cres.setdefault("sex_effect", {})[feat] = entry

    strat_stats[name] = dict(sex_source=c["sex_source"], note=c["note"], results=cres)

    # ---- interaction: expression ~ group * sex ----
    ires = {}
    if c["ctrl_label"] is not None:
        nF = (fdf["sex"] == "Female").sum()
        nM = (fdf["sex"] == "Male").sum()
        if nF >= MIN_INTERACTION_N and nM >= MIN_INTERACTION_N:
            d = fdf.copy()
            d["group"] = (d["group"] == c["case_label"]).astype(int)
            d["sexM"] = (d["sex"] == "Male").astype(int)
            for feat in fdf.columns:
                if feat in ["group", "sex"]:
                    continue
                dd = d[[feat, "group", "sexM"]].dropna()
                try:
                    fit = smf.ols(f"Q('{feat}') ~ group * sexM", data=dd).fit()
                    key = [k for k in fit.pvalues.index if ":" in k][0]
                    ires[feat] = dict(
                        interaction_coef=float(fit.params[key]),
                        interaction_p=float(fit.pvalues[key]),
                        group_p=float(fit.pvalues.get("group", np.nan)),
                        sex_p=float(fit.pvalues.get("sexM", np.nan)),
                        n=len(dd), n_F=int(nF), n_M=int(nM),
                    )
                except Exception as ex:
                    ires[feat] = dict(error=str(ex))
        else:
            ires["skipped"] = f"sex stratum <{MIN_INTERACTION_N} (F={nF}, M={nM})"
    interaction_tests[name] = ires

    # collect plot data (axis scores only)
    for feat in ["LDAM", "Ferroptosis", "WGCNA_brown_module"]:
        if feat in fdf.columns:
            for _, row in fdf[[feat, "group", "sex"]].iterrows():
                plot_data.append(dict(cohort=name, score=feat, value=row[feat],
                                      group=row["group"], sex=row["sex"]))

with open(OUT / "stratified_stats.json", "w") as f:
    json.dump(strat_stats, f, indent=2)
with open(OUT / "interaction_tests.json", "w") as f:
    json.dump(interaction_tests, f, indent=2)
print("stats written")

# ---------------- figure ----------------
pd_plot = pd.DataFrame(plot_data)
score_labels = {"LDAM": "LDAM score", "Ferroptosis": "Ferroptosis score",
                "WGCNA_brown_module": "WGCNA brown module score"}
fig_cohorts = [n for n in ["GSE134697", "GSE168375", "GSE143272"]
               if n in pd_plot["cohort"].unique()]
scores = ["LDAM", "Ferroptosis", "WGCNA_brown_module"]
fig, axes = plt.subplots(len(scores), len(fig_cohorts),
                         figsize=(3.4 * len(fig_cohorts), 3.2 * len(scores)), squeeze=False)
sex_colors = {"Female": "#C44E52", "Male": "#4C72B0"}
offsets = {"Female": -0.18, "Male": 0.18}
rng = np.random.default_rng(42)

for j, coh in enumerate(fig_cohorts):
    for i, sc in enumerate(scores):
        ax = axes[i][j]
        sub = pd_plot[(pd_plot["cohort"] == coh) & (pd_plot["score"] == sc)]
        groups = [g for g in ["Control", "TLE"] if g in sub["group"].unique()]
        for gi, grp in enumerate(groups):
            for sexv in ["Female", "Male"]:
                vals = sub[(sub["group"] == grp) & (sub["sex"] == sexv)]["value"].values
                if len(vals) == 0:
                    continue
                pos = gi + offsets[sexv]
                bp = ax.boxplot(vals, positions=[pos], widths=0.3, showfliers=False,
                                patch_artist=True, medianprops=dict(color="black", lw=1.2),
                                whiskerprops=dict(color=sex_colors[sexv]),
                                capprops=dict(color=sex_colors[sexv]))
                bp["boxes"][0].set_facecolor(sex_colors[sexv])
                bp["boxes"][0].set_alpha(0.45)
                bp["boxes"][0].set_edgecolor(sex_colors[sexv])
                ax.scatter(pos + rng.uniform(-0.08, 0.08, len(vals)), vals,
                           s=9, color=sex_colors[sexv], alpha=0.75, zorder=3,
                           edgecolors="none")
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, fontsize=9)
        if i == 0:
            src = cohorts[coh]["sex_source"]
            ax.set_title(f"{coh}\n(sex: {src})", fontsize=9)
        if j == 0:
            ax.set_ylabel(score_labels[sc], fontsize=9)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

handles = [plt.Line2D([], [], marker="s", ls="", color=sex_colors[s], label=s)
           for s in ["Female", "Male"]]
fig.legend(handles=handles, loc="upper right", fontsize=9, frameon=False)
fig.suptitle("Sex-stratified axis scores: TLE vs control across bulk cohorts", fontsize=12)
fig.tight_layout(rect=[0, 0, 0.95, 0.96])
fig.savefig(FIG / "Fig_sex_stratified.png", dpi=300)
print("figure written:", FIG / "Fig_sex_stratified.png")
