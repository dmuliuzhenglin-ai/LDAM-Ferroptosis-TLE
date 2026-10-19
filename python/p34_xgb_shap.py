# -*- coding: utf-8 -*-
"""
p34_xgb_shap.py
XGBoost retraining + SHAP interpretability of the 7-gene diagnostic signature
(LPCAT2, S100A10, CD9, SYNGR2 | ITGAM, ITGB2, YBX1) on GSE143272 blood
(34 drug-naive TLE vs 51 healthy controls), external validation on GSE168375.

Provenance:
  - Original 4-gene signature: LASSO+SVM-RFE+RF intersection
    (results/04_ML_selected_genes.csv, R/04_machine_learning_hub_ROC.R)
  - Blood exosome hub genes: results/33_blood_exosome/33_blood_hub_genes.csv
  - Combined 7-gene logistic model: results/37_blood_combined_ROC/
    (apparent AUC=0.803 [0.706-0.900], stratified 10-fold x10 CV AUC=0.759+/-0.023)
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, confusion_matrix
import xgboost as xgb
import shap

BASE = Path(r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究")
OUT = BASE / "results" / "69_xgb_shap"
FIG = BASE / "figures" / "v25"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GENES7 = ["LPCAT2", "S100A10", "CD9", "SYNGR2", "ITGAM", "ITGB2", "YBX1"]
SIG_ORIG = ["LPCAT2", "S100A10", "CD9", "SYNGR2"]
SIG_BLOOD = ["ITGAM", "ITGB2", "CD9", "YBX1"]

XGB_PARAMS = dict(
    n_estimators=200, max_depth=2, learning_rate=0.05,
    subsample=0.8, colsample_bytree=0.8, reg_lambda=1.0,
    eval_metric="logloss", random_state=0, n_jobs=4,
)


def load_cohort(expr_file, meta_file):
    e = pd.read_csv(BASE / expr_file)
    e = e.rename(columns={e.columns[0]: "gene"})
    e = e.groupby("gene", as_index=True).mean(numeric_only=True)
    m = pd.read_csv(BASE / meta_file)
    samples = [s for s in m["sample"] if s in e.columns]
    grp = m.set_index("sample").loc[samples, "group"]
    y = (grp.str.lower() == "epilepsy").astype(int).values
    mat = e.loc[[g for g in GENES7 if g in e.index], samples]
    # z-score per gene across samples (within cohort), as in analysis 37
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0)
    z = z.fillna(0.0)
    X = z.T[GENES7]  # samples x 7, fixed column order
    return X, y, samples


X, y, samples = load_cohort(
    "results/08_GSE143272_normalized.csv",
    "results/08_GSE143272_sample_meta.csv")
print(f"GSE143272: {X.shape[0]} samples, {y.sum()} TLE vs {(y == 0).sum()} control")

Xv, yv, samples_v = load_cohort(
    "results/02_GSE168375_vst_normalized.csv",
    "results/02_GSE168375_sample_meta.csv")
print(f"GSE168375 (external): {Xv.shape[0]} samples, {yv.sum()} epilepsy vs {(yv == 0).sum()} control")

# ------------------------------------------------------------------
# 1. Repeated 5-fold CV (10 different seeds => 50 fits)
# ------------------------------------------------------------------
N_SEEDS = 10
cv_aucs, cv_accs, cv_sens, cv_spec = [], [], [], []
oof_prob_sum = np.zeros(len(y))
oof_prob_cnt = np.zeros(len(y))
rng_seeds = [2024 + i for i in range(N_SEEDS)]

for seed in rng_seeds:
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for tr, te in skf.split(X, y):
        clf = xgb.XGBClassifier(**XGB_PARAMS)
        clf.fit(X.iloc[tr], y[tr])
        p = clf.predict_proba(X.iloc[te])[:, 1]
        cv_aucs.append(roc_auc_score(y[te], p))
        pred = (p >= 0.5).astype(int)
        cv_accs.append(accuracy_score(y[te], pred))
        tn, fp, fn, tp = confusion_matrix(y[te], pred, labels=[0, 1]).ravel()
        cv_sens.append(tp / (tp + fn))
        cv_spec.append(tn / (tn + fp))
        oof_prob_sum[te] += p
        oof_prob_cnt[te] += 1

oof_prob = oof_prob_sum / oof_prob_cnt
oof_auc = roc_auc_score(y, oof_prob)
oof_pred = (oof_prob >= 0.5).astype(int)
tn, fp, fn, tp = confusion_matrix(y, oof_pred, labels=[0, 1]).ravel()

cv_perf = {
    "design": "5-fold stratified CV x 10 repeats (seeds 2024-2033), 50 fits",
    "n_samples": int(len(y)), "n_TLE": int(y.sum()), "n_control": int((y == 0).sum()),
    "cv_auc_mean": float(np.mean(cv_aucs)), "cv_auc_sd": float(np.std(cv_aucs, ddof=1)),
    "cv_accuracy_mean": float(np.mean(cv_accs)), "cv_accuracy_sd": float(np.std(cv_accs, ddof=1)),
    "cv_sensitivity_mean": float(np.mean(cv_sens)), "cv_sensitivity_sd": float(np.std(cv_sens, ddof=1)),
    "cv_specificity_mean": float(np.mean(cv_spec)), "cv_specificity_sd": float(np.std(cv_spec, ddof=1)),
    "oof_pooled_auc": float(oof_auc),
    "oof_accuracy": float((oof_pred == y).mean()),
    "oof_sensitivity": float(tp / (tp + fn)), "oof_specificity": float(tn / (tn + fp)),
    "cv_aucs_all": [float(a) for a in cv_aucs],
}

# ------------------------------------------------------------------
# 2. Full-data model + apparent performance + external validation
# ------------------------------------------------------------------
final_clf = xgb.XGBClassifier(**XGB_PARAMS)
final_clf.fit(X, y)
app_prob = final_clf.predict_proba(X)[:, 1]
app_auc = roc_auc_score(y, app_prob)
val_prob = final_clf.predict_proba(Xv)[:, 1]
val_auc = roc_auc_score(yv, val_prob)
val_pred = (val_prob >= 0.5).astype(int)
tnv, fpv, fnv, tpv = confusion_matrix(yv, val_pred, labels=[0, 1]).ravel()

cv_perf["apparent_auc_train"] = float(app_auc)
cv_perf["external_validation"] = {
    "cohort": "GSE168375 (hippocampus tissue, epilepsy vs control)",
    "n": int(len(yv)), "n_epilepsy": int(yv.sum()), "n_control": int((yv == 0).sum()),
    "auc": float(val_auc),
    "sensitivity": float(tpv / (tpv + fnv)), "specificity": float(tnv / (tnv + fpv)),
}
cv_perf["prior_benchmark_analysis37"] = {
    "model": "7-gene logistic regression (glm), same cohort GSE143272",
    "apparent_auc": 0.803344867358708, "apparent_ci": [0.70636897262635, 0.900320762091066],
    "cv_auc": 0.759, "cv_sd": 0.0231400726654642,
    "cv_design": "stratified 10-fold x 10 repeats",
    "source": "results/37_blood_combined_ROC/37_roc_models.csv",
}

# ------------------------------------------------------------------
# 3. SHAP
# ------------------------------------------------------------------
# shap 0.49.1 cannot parse xgboost 3.x vector base_score -> use the booster's
# native pred_contribs / pred_interactions (identical TreeSHAP output).
booster = final_clf.get_booster()
dm = xgb.DMatrix(X, feature_names=GENES7)
contrib = booster.predict(dm, pred_contribs=True)      # (n, 7+1) last col = bias
shap_vals = contrib[:, :-1]                            # (n, 7)
inter_full = booster.predict(dm, pred_interactions=True)  # (n, 8, 8)
inter_vals = inter_full[:, :-1, :-1]                   # (n, 7, 7)

mean_abs = np.abs(shap_vals).mean(axis=0)
order = np.argsort(-mean_abs)
ranking = [GENES7[i] for i in order]

# direction: correlation between expression (z) and SHAP value per gene
deg = pd.read_csv(BASE / "results/08_GSE143272_DEG_Epilepsy_vs_Control.csv")
deg = deg.set_index("gene")
direction_check = {}
for gi, g in enumerate(GENES7):
    r = np.corrcoef(X[g].values, shap_vals[:, gi])[0, 1]
    shap_dir = "up_pushes_TLE" if r > 0 else "up_pushes_control"
    lfc = float(deg.loc[g, "log2FoldChange"]) if g in deg.index else None
    deg_dir = "up_in_TLE" if (lfc is not None and lfc > 0) else "down_in_TLE"
    consistent = (r > 0) == (lfc > 0)
    direction_check[g] = {
        "corr_expr_shap": float(r),
        "shap_direction": shap_dir,
        "deg_log2FC_epilepsy_vs_control": lfc,
        "deg_direction": deg_dir,
        "consistent": bool(consistent),
        "mean_abs_shap": float(mean_abs[gi]),
        "rank_by_mean_abs_shap": int(ranking.index(g) + 1),
    }

# strongest interaction (off-diagonal mean |interaction|)
offdiag = np.abs(inter_vals).mean(axis=0)
np.fill_diagonal(offdiag, 0.0)
i1, i2 = np.unravel_index(np.argmax(offdiag), offdiag.shape)
top_pair = (GENES7[i1], GENES7[i2])
top_pair_strength = float(offdiag[i1, i2])

# ------------------------------------------------------------------
# 4. Save results
# ------------------------------------------------------------------
with open(OUT / "cv_performance.json", "w", encoding="utf-8") as f:
    json.dump(cv_perf, f, indent=2, ensure_ascii=False)

sv = pd.DataFrame(shap_vals, columns=[f"shap_{g}" for g in GENES7])
sv.insert(0, "sample", samples)
sv["label"] = y
for g in GENES7:
    sv[f"expr_z_{g}"] = X[g].values
sv.to_parquet(OUT / "shap_values.parquet", index=False)
np.savez_compressed(OUT / "shap_interactions.npz",
                    interactions=inter_vals, genes=np.array(GENES7),
                    samples=np.array(samples))

gene_direction = {
    "cohort": "GSE143272 blood, drug-naive TLE vs healthy",
    "deg_source": "results/08_GSE143272_DEG_Epilepsy_vs_Control.csv",
    "genes": direction_check,
    "strongest_interaction_pair": list(top_pair),
    "interaction_strength_mean_abs": top_pair_strength,
}
with open(OUT / "gene_direction_check.json", "w", encoding="utf-8") as f:
    json.dump(gene_direction, f, indent=2, ensure_ascii=False)

# ------------------------------------------------------------------
# 5. Figures
# ------------------------------------------------------------------
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False})

# --- Fig 1: SHAP summary (beeswarm + bar + dependence) ---
fig = plt.figure(figsize=(14, 5.0), dpi=300)
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 0.95, 1.1], wspace=0.5)

fig.add_subplot(gs[0])
shap.summary_plot(shap_vals, X, feature_names=GENES7, show=False, plot_size=None)
ax1 = plt.gca()
ax1.set_xlabel("SHAP value")
ax1.set_title("A  SHAP beeswarm (GSE143272)", loc="left", fontweight="bold")

fig.add_subplot(gs[1])
shap.summary_plot(shap_vals, X, feature_names=GENES7, plot_type="bar",
                  show=False, plot_size=None, color="#1783FE")
ax2 = plt.gca()
ax2.set_xlabel("Mean |SHAP value|")
ax2.set_title("B  Mean |SHAP| importance", loc="left", fontweight="bold")

ax3 = fig.add_subplot(gs[2])
g1, g2 = top_pair
sc = ax3.scatter(X[g1], shap_vals[:, GENES7.index(g1)],
                 c=X[g2], cmap="coolwarm", s=22, edgecolor="grey", linewidth=0.3)
ax3.axhline(0, color="grey", lw=0.6, ls="--")
ax3.set_xlabel(f"{g1} expression (z-score)")
ax3.set_ylabel(f"SHAP value for {g1}")
cb = fig.colorbar(sc, ax=ax3)
cb.set_label(f"{g2} expression (z)")
ax3.set_title(f"C  Dependence: {g1} x {g2}\n(strongest interaction)",
              loc="left", fontweight="bold")
fig.savefig(FIG / "Fig_shap_summary.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# --- Fig 2: performance (CV distribution + ROC + comparison) ---
fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), dpi=300)

ax = axes[0]
ax.boxplot(cv_aucs, positions=[0], widths=0.35, patch_artist=True,
           boxprops=dict(facecolor="#DCE9FD"), medianprops=dict(color="#D64545"))
ax.scatter(np.zeros(len(cv_aucs)) + np.random.uniform(-0.06, 0.06, len(cv_aucs)),
           cv_aucs, s=10, color="#1783FE", alpha=0.6, zorder=3)
prior_cv = cv_perf["prior_benchmark_analysis37"]["cv_auc"]
ax.axhline(prior_cv, color="#67B8A6", ls="--", lw=1.5,
           label=f"Prior 7-gene logit CV AUC = {prior_cv:.3f}")
ax.set_xticks([0]); ax.set_xticklabels(["XGBoost"])
ax.set_ylabel("CV AUC (50 folds: 5-fold x 10 seeds)")
ax.set_title("A  Repeated CV AUC distribution", loc="left", fontweight="bold")
ax.legend(loc="lower right", fontsize=8)

ax = axes[1]
fpr, tpr, _ = roc_curve(y, oof_prob)
ax.plot(fpr, tpr, color="#1783FE", lw=2,
        label=f"OOF pooled (AUC = {oof_auc:.3f})")
fv, tv, _ = roc_curve(yv, val_prob)
ax.plot(fv, tv, color="#D64545", lw=2,
        label=f"External GSE168375 (AUC = {val_auc:.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=0.8)
ax.set_xlabel("1 - Specificity"); ax.set_ylabel("Sensitivity")
ax.set_title("B  ROC curves", loc="left", fontweight="bold")
ax.legend(loc="lower right", fontsize=8)

ax = axes[2]
labels = ["XGBoost\nCV (5Fx10)", "Prior logit\nCV (10Fx10)",
          "XGBoost\napparent", "Prior logit\napparent", "XGBoost\nexternal\nGSE168375"]
vals = [np.mean(cv_aucs), prior_cv, app_auc,
        cv_perf["prior_benchmark_analysis37"]["apparent_auc"], val_auc]
errs = [np.std(cv_aucs, ddof=1), cv_perf["prior_benchmark_analysis37"]["cv_sd"],
        0, 0, 0]
colors = ["#1783FE", "#67B8A6", "#8C98A8", "#B9C2CE", "#D64545"]
ax.bar(range(len(vals)), vals, yerr=errs, capsize=4, color=colors)
for i, v in enumerate(vals):
    ax.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=8)
ax.set_xticks(range(len(vals))); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylim(0, 1.05); ax.set_ylabel("AUC")
ax.set_title("C  XGBoost vs prior logistic benchmark", loc="left", fontweight="bold")

fig.tight_layout()
fig.savefig(FIG / "Fig_xgb_performance.png", bbox_inches="tight", dpi=300)
plt.close(fig)

print(json.dumps({"cv_auc_mean": cv_perf["cv_auc_mean"],
                  "cv_auc_sd": cv_perf["cv_auc_sd"],
                  "oof_auc": oof_auc, "apparent_auc": app_auc,
                  "external_auc": val_auc,
                  "ranking": ranking,
                  "top_pair": top_pair,
                  "pair_strength": top_pair_strength}, indent=2))
print("DONE -> results/69_xgb_shap/ and figures/v25/")
