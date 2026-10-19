"""A3: pySCENIC regulon analysis on GSE190452 microglia.
GRNBoost2 (arboreto) -> cisTarget pruning (mc9nr) -> AUCell ->
differential regulon activity: LDAM-high vs LDAM-low microglia (quartiles) and HTLE vs Control.
Run with gsmap_env python. Stages: grn, ctx, aucell, diff (default: all).
"""
import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
OUT = os.path.join(BASE, "results", "56_pySCENIC")
os.makedirs(OUT, exist_ok=True)
CIS = os.path.join(BASE, "data", "cistarget")
FEATHERS = [os.path.join(CIS, "hg38_500bp.mc9nr.genes_vs_motifs.rankings.feather"),
            os.path.join(CIS, "hg38_10kb.mc9nr.genes_vs_motifs.rankings.feather")]
MOTIF2TF = os.path.join(CIS, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
TF_LIST = os.path.join(CIS, "allTFs_hg38.txt")

stage = sys.argv[1] if len(sys.argv) > 1 else "all"

def main():
    _run(stage)

def _run(stage):
    # ---- expression prep (cached) ---
    EX_MTX = os.path.join(OUT, "microglia_counts_filtered.tsv.gz")
    if stage in ("grn", "all"):
        a = sc.read_h5ad(os.path.join(BASE, "results", "43_scDRS", "gse190452_microglia.h5ad"))
        counts = a.layers["counts"]
        if sp.issparse(counts):
            counts = counts.tocsr()
        # gene filter: expressed (>=1 count) in >=1% of cells
        frac = np.asarray((counts > 0).sum(axis=0)).ravel() / counts.shape[0]
        keep_gene = frac >= 0.01
        tfs = set(pd.read_csv(TF_LIST, header=None)[0].str.strip())
        keep_tf = np.array([g in tfs for g in a.var_names])
        keep = keep_gene | keep_tf
        ex = pd.DataFrame(counts[:, keep].toarray() if sp.issparse(counts[:, keep]) else counts[:, keep],
                          index=a.obs_names, columns=a.var_names[keep])
        ex.to_csv(EX_MTX, sep="\t", compression="gzip")
        print("expr matrix:", ex.shape, flush=True)
        del a, counts, ex

    if stage in ("grn", "all"):
        from arboreto.algo import grnboost2
        from distributed import LocalCluster, Client
        ex = pd.read_csv(EX_MTX, sep="\t", index_col=0)
        tf_names = [t for t in pd.read_csv(TF_LIST, header=None)[0].str.strip() if t in ex.columns]
        cluster = LocalCluster(n_workers=1, threads_per_worker=8, memory_limit="24GB",
                               processes=False, silence_logs=40)
        client = Client(cluster)
        adj = grnboost2(expression_data=ex, tf_names=tf_names, client_or_address=client, seed=0)
        client.close(); cluster.close()
        adj.to_csv(os.path.join(OUT, "grn_adjacencies.tsv"), sep="\t", index=False, header=False)
        print("GRN edges:", len(adj), flush=True)

    if stage in ("ctx", "all"):
        from pyscenic.utils import modules_from_adjacencies
        from pyscenic.prune import prune2df, df2regulons
        ex = pd.read_csv(EX_MTX, sep="\t", index_col=0)
        adj = pd.read_csv(os.path.join(OUT, "grn_adjacencies.tsv"), sep="\t", header=None,
                          names=["TF", "target", "importance"])
        modules = modules_from_adjacencies(adj, ex)
        from distributed import LocalCluster, Client
        cluster = LocalCluster(n_workers=4, threads_per_worker=1, memory_limit="5GB",
                               processes=True, silence_logs=40)
        client = Client(cluster)
        from ctxcore.rnkdb import FeatherRankingDatabase
        rnkdbs = [FeatherRankingDatabase(fname=f, name=os.path.basename(f)) for f in FEATHERS]
        df = prune2df(rnkdbs, modules, MOTIF2TF, client_or_address=client)
        client.close(); cluster.close()
        df.to_csv(os.path.join(OUT, "regulons_df.tsv"), sep="\t")
        # build regulons manually: df is indexed by TF, columns MultiIndex ("Enrichment", field);
        # TargetGenes cells may be strings in this pySCENIC/dask version -> literal_eval
        import ast, pickle
        from ctxcore.genesig import Regulon
        d = df.copy()
        if d.columns.nlevels == 2:
            d.columns = d.columns.droplevel(0)
        regs = []
        for tf_name, grp in d.groupby(level=0):
            targets = {}
            activating = True
            for _, r in grp.iterrows():
                ctx_str = str(r["Context"])
                if "repressing" in ctx_str:
                    activating = False
                tg = r["TargetGenes"]
                if isinstance(tg, str):
                    tg = ast.literal_eval(tg)
                for gene, w in tg:
                    if gene not in targets or w > targets[gene]:
                        targets[gene] = w
            if len(targets) >= 5:
                regs.append(Regulon(name=f"{tf_name}({'+' if activating else '-'})",
                                    gene2weight=targets))
        with open(os.path.join(OUT, "regulons.pkl"), "wb") as f:
            pickle.dump(regs, f)
        print("regulons:", len(regs), flush=True)

    if stage in ("aucell", "all"):
        from pyscenic.aucell import aucell
        import pickle
        with open(os.path.join(OUT, "regulons.pkl"), "rb") as f:
            regs = pickle.load(f)
        ex = pd.read_csv(EX_MTX, sep="\t", index_col=0)
        auc = aucell(ex, regs, num_workers=8, seed=0)  # cells x regulons
        auc.to_csv(os.path.join(OUT, "aucell_matrix.csv.gz"), compression="gzip")
        print("AUCell:", auc.shape, flush=True)

    if stage in ("diff", "all"):
        auc = pd.read_csv(os.path.join(OUT, "aucell_matrix.csv.gz"), index_col=0)
        ann = pd.read_csv(os.path.join(BASE, "results", "43_scDRS", "cell_annotations.csv"), index_col=0)
        ann = ann.loc[ann["celltype"] == "Microglia"]
        auc = auc.loc[auc.index.intersection(ann.index)]
        ann = ann.loc[auc.index]
        auc.columns = [c.replace("(", "").replace(")", "").replace(" ", "_") for c in auc.columns]
        from scipy import stats as st
        rows = []
        q75, q25 = ann["LDAM_score"].quantile(0.75), ann["LDAM_score"].quantile(0.25)
        hi = ann["LDAM_score"] >= q75
        lo = ann["LDAM_score"] <= q25
        cond = ann["condition"] == "HTLE"
        for reg in auc.columns:
            v = auc[reg]
            a, b = v[hi], v[lo]
            t1, p1 = st.mannwhitneyu(a, b, alternative="greater")[:2] if False else (np.nan,)*2
            u, p1 = st.mannwhitneyu(a, b, alternative="greater")
            rbc1 = 1 - 2 * u / (len(a) * len(b))  # rank-biserial (a>b positive)
            h_, c_ = v[cond], v[~cond]
            u2, p2 = st.mannwhitneyu(h_, c_, alternative="greater")
            rbc2 = 1 - 2 * u2 / (len(h_) * len(c_))
            rows.append(dict(regulon=reg, n_hi=len(a), n_lo=len(b),
                             auc_hi=a.mean(), auc_lo=b.mean(), delta_ldam=a.mean() - b.mean(),
                             p_ldam=p1, rbc_ldam=rbc1,
                             auc_htle=h_.mean(), auc_ctrl=c_.mean(), delta_htle=h_.mean() - c_.mean(),
                             p_htle=p2, rbc_htle=rbc2))
        res = pd.DataFrame(rows)
        from statsmodels.stats.multitest import multipletests
        res["fdr_ldam"] = multipletests(res["p_ldam"], method="fdr_bh")[1]
        res["fdr_htle"] = multipletests(res["p_htle"], method="fdr_bh")[1]
        res = res.sort_values("p_ldam")
        res.to_csv(os.path.join(OUT, "regulon_differential.csv"), index=False)
        print(res.head(20).to_string(index=False), flush=True)
    print("PYSCENIC_STAGE_DONE", stage, flush=True)


if __name__ == "__main__":
    main()
