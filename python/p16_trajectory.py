"""Task 4: microglia state trajectory — diffusion pseudotime + LDAM/ferroptosis coupling.
Data: results/43_scDRS/gse190452_microglia.h5ad (5982 microglia, layers counts).
Outputs: results/58_trajectory/
"""
import os, sys
import numpy as np, pandas as pd
import scanpy as sc
from scipy import stats as st

BASE = os.getcwd()
OUT = os.path.join(BASE, "results", "58_trajectory")
os.makedirs(OUT, exist_ok=True)

def main():
    ad = sc.read_h5ad(os.path.join(BASE, "results", "43_scDRS", "gse190452_microglia.h5ad"))
    ann = pd.read_csv(os.path.join(BASE, "results", "43_scDRS", "cell_annotations.csv"), index_col=0)
    ann = ann.loc[ann["celltype"] == "Microglia"]
    ad = ad[ad.obs_names.intersection(ann.index)].copy()
    ann = ann.loc[ad.obs_names]
    ad.obs["condition"] = ann["condition"]
    ad.obs["LDAM_score"] = ann["LDAM_score"]
    ad.obs["Ferroptosis_score"] = ann["Ferroptosis_score"]

    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=3000, flavor="seurat")
    sc.pp.scale(ad, max_value=10)
    sc.pp.pca(ad, n_comps=40)
    sc.pp.neighbors(ad, n_neighbors=20, n_pcs=30)
    sc.tl.diffmap(ad)
    # root: cell with lowest LDAM score
    root = int(np.argmin(ad.obs["LDAM_score"].values))
    ad.uns["iroot"] = root
    sc.tl.dpt(ad)
    sc.tl.umap(ad)
    ad.obs["dpt"] = ad.obs["dpt_pseudotime"]

    # correlations with scores
    res = {}
    for s in ["LDAM_score", "Ferroptosis_score"]:
        rho, p = st.spearmanr(ad.obs["dpt"], ad.obs[s])
        res[s] = (rho, p)
    # condition distribution along pseudotime (terciles)
    t = pd.qcut(ad.obs["dpt"], 3, labels=["early", "mid", "late"])
    tab = pd.crosstab(t, ad.obs["condition"], normalize="index")
    chi2, pchi, _, _ = st.chi2_contingency(pd.crosstab(t, ad.obs["condition"]))
    # gene trends along dpt (top genes correlated with dpt)
    hv = ad.var_names[ad.var["highly_variable"]]
    X = ad[:, hv].X
    X = X.toarray() if hasattr(X, "toarray") else X
    dpt = ad.obs["dpt"].values
    cors = np.array([st.spearmanr(dpt, X[:, i])[0] for i in range(X.shape[1])])
    gene_tr = pd.DataFrame({"gene": hv, "rho_dpt": cors}).dropna().sort_values("rho_dpt")
    gene_tr["direction"] = np.where(gene_tr.rho_dpt > 0, "late_up", "early_up")

    # LDAM / ferroptosis genes among late-up
    ldam = set(pd.read_csv(os.path.join(BASE, "results", "03_WGCNA_LDAM_hub_module_genes.txt"), header=None)[0])
    fer = pd.read_csv(os.path.join(BASE, "results", "05_KEGG_candidates.csv"))
    late_up = set(gene_tr[gene_tr.rho_dpt > 0.15].gene)
    n_late = len(late_up); n_all = len(hv)
    ov_ldam = len(late_up & ldam)
    from scipy.stats import fisher_exact
    _, p_ldam = fisher_exact([[ov_ldam, n_late - ov_ldam], [len(ldam & set(hv)) - ov_ldam, n_all - n_late - (len(ldam & set(hv)) - ov_ldam)]])

    ad.write_h5ad(os.path.join(OUT, "microglia_dpt.h5ad"))
    gene_tr.to_csv(os.path.join(OUT, "gene_pseudotime_trends.csv"), index=False)
    summ = {
        "n_cells": ad.n_obs,
        "spearman_dpt_LDAM": res["LDAM_score"],
        "spearman_dpt_ferroptosis": res["Ferroptosis_score"],
        "condition_x_tercile": tab.to_dict(),
        "chi2_p": pchi,
        "late_up_genes": n_late,
        "late_up_in_LDAM_set": ov_ldam,
        "fisher_p_LDAM_late": p_ldam,
        "top_late": gene_tr.tail(15).gene.tolist(),
        "top_early": gene_tr.head(15).gene.tolist(),
    }
    import json
    json.dump(summ, open(os.path.join(OUT, "trajectory_summary.json"), "w"), indent=2, default=str)
    print(json.dumps(summ, indent=2, default=str))

if __name__ == "__main__":
    main()
