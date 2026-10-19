"""Task 5: metacell co-expression modules in microglia (hdWGCNA-style, Python).
Metacells by kNN aggregation -> Pearson corr on HVGs -> hierarchical modules.
Module eigengenes tested against LDAM/ferroptosis scores, condition, and IKZF1/IRF8 AUCell.
Outputs: results/59_metacell_wgcna/
"""
import os
import numpy as np, pandas as pd
import scanpy as sc
from scipy import stats as st
from scipy.cluster import hierarchy as hc

BASE = os.getcwd()
OUT = os.path.join(BASE, "results", "59_metacell_wgcna")
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

    X = ad.layers["counts"]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    lib = X.sum(1, keepdims=True)
    Xn = np.log1p(X / np.maximum(lib, 1) * 1e4)

    # HVGs by variance of log-normalized
    v = Xn.var(0)
    hvg_idx = np.argsort(-v)[:2000]
    genes = np.array(ad.var_names)[hvg_idx]
    Xh = Xn[:, hvg_idx]

    # metacells: kNN on PCA, aggregate 25 nearest cells (excluding self)
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    Z = PCA(n_components=30, random_state=0).fit_transform((Xh - Xh.mean(0)) / (Xh.std(0) + 1e-8))
    k = 25
    nn = NearestNeighbors(n_neighbors=k).fit(Z)
    _, idx = nn.kneighbors(Z)
    rng = np.random.default_rng(0)
    # sample 150 metacells from random seeds without heavy overlap control
    seeds = rng.choice(len(Z), size=150, replace=False)
    meta = np.stack([Xh[idx[s]].mean(0) for s in seeds])
    meta_cond = np.array([ann["condition"].values[s] for s in seeds])
    meta_ldam = np.array([ann["LDAM_score"].values[idx[s]].mean() for s in seeds])
    meta_fer = np.array([ann["Ferroptosis_score"].values[idx[s]].mean() for s in seeds])

    # gene-gene correlation across metacells -> modules via hierarchical clustering
    C = np.corrcoef(meta.T)
    C = np.nan_to_num(C)
    dist = 1 - np.clip(C, -1, 1)
    np.fill_diagonal(dist, 0)
    dist = (dist + dist.T) / 2
    link = hc.linkage(dist[np.triu_indices(len(genes), 1)], method="average")
    labels = hc.fcluster(link, t=12, criterion="maxclust")
    mod = pd.DataFrame({"gene": genes, "module": labels})
    counts = mod.module.value_counts()
    keep_mods = counts[counts >= 20].index
    mod = mod[mod.module.isin(keep_mods)]
    # eigengene = PC1 of module expression across metacells
    rows = []
    eig = {}
    for m in sorted(keep_mods):
        g = mod[mod.module == m].gene.values
        gi = [np.where(genes == x)[0][0] for x in g]
        M = meta[:, gi]
        Mz = (M - M.mean(0)) / (M.std(0) + 1e-8)
        u, s_, vt = np.linalg.svd(Mz, full_matrices=False)
        e1 = u[:, 0] * s_[0]
        # orient: positive correlation with module mean
        if np.corrcoef(e1, Mz.mean(1))[0, 1] < 0:
            e1 = -e1
        eig[m] = e1
        for name, vec in [("LDAM", meta_ldam), ("Ferroptosis", meta_fer)]:
            rho, p = st.spearmanr(e1, vec)
            rows.append(dict(module=int(m), n_genes=len(g), trait=name, rho=rho, p=p))
        # condition difference of eigengene
        a = e1[meta_cond == "HTLE"]; b = e1[meta_cond == "Control"]
        if len(a) > 2 and len(b) > 2:
            u_, p_ = st.mannwhitneyu(a, b, alternative="two-sided")
            rows.append(dict(module=int(m), n_genes=len(g), trait="HTLE_vs_Control",
                             rho=float(a.mean() - b.mean()), p=p_))
    df = pd.DataFrame(rows)
    df["fdr"] = df.groupby("trait").p.transform(lambda p: st.false_discovery_control(p) if hasattr(st, "false_discovery_control") else p)
    mod.to_csv(os.path.join(OUT, "metacell_modules.csv"), index=False)
    df.to_csv(os.path.join(OUT, "module_trait_associations.csv"), index=False)
    # AUCell overlay: mean regulon activity per metacell for top regulons
    auc = pd.read_csv(os.path.join(BASE, "results", "56_pySCENIC", "aucell_matrix.csv.gz"), index_col=0)
    auc = auc.loc[ad.obs_names]
    reg_rows = []
    for reg in ["IKZF1(+)", "ETV6(+)", "JAZF1(+)", "RUNX1(+)", "IRF8(+)", "FLI1(+)", "FOXP2(+)", "ETV5(+)"]:
        if reg not in auc.columns: continue
        v = auc[reg].values
        mv = np.array([v[idx[s]].mean() for s in seeds])
        for m in sorted(keep_mods):
            rho, p = st.spearmanr(eig[m], mv)
            reg_rows.append(dict(module=int(m), regulon=reg, rho=rho, p=p))
    pd.DataFrame(reg_rows).to_csv(os.path.join(OUT, "module_regulon_overlay.csv"), index=False)
    print(df.sort_values(["trait", "p"]).groupby("trait").head(5).to_string(index=False))
    print("modules kept:", sorted(keep_mods), "metacells:", len(seeds))

if __name__ == "__main__":
    main()
