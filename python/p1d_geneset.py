"""Phase 1 step D: competitive gene-set analysis (MAGMA-style) + hub gene table."""
import numpy as np, pandas as pd, os, glob
from scipy.stats import norm, t as tdist

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
RESD = os.path.join(BASE, "results", "41_MAGMA")

gp = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(RESD, "gene_p_chr*.csv"))], ignore_index=True)
gp["P"] = gp["P"].clip(1e-300, 1.0)
gp["Z"] = norm.isf(gp["P"])
gp = gp[np.isfinite(gp["Z"])].reset_index(drop=True)
print("total genes tested:", len(gp))

sets = {}
sets["LDAM"] = pd.read_csv(os.path.join(BASE, "genesets", "LDAM_geneset.txt"), header=None)[0].str.strip().tolist()
sets["Ferroptosis"] = pd.read_csv(os.path.join(BASE, "genesets", "Ferroptosis_geneset.txt"), header=None)[0].str.strip().tolist()
bm = pd.read_csv(os.path.join(BASE, "results", "03_WGCNA_LDAM_hub_module_genes.csv"))
sets["WGCNA_brown"] = bm["gene"].astype(str).str.strip().tolist()

rows = []
for name, members in sets.items():
    ind = gp["symbol"].isin(set(members)).astype(float).values
    n_in = int(ind.sum())
    X = np.column_stack([np.ones(len(gp)), ind, np.log(gp["nSNP"].values)])
    y = gp["Z"].values
    beta, res_, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    s2 = resid @ resid / dof
    covb = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(covb[1, 1])
    tstat = beta[1] / se
    p_one = tdist.sf(tstat, dof)  # direction=greater
    rows.append(dict(set=name, n_genes_in_set=len(members), n_tested=n_in,
                     BETA=beta[1], SE=se, T=tstat, P=p_one))
    print(name, rows[-1])

res = pd.DataFrame(rows).sort_values("P")
p = res["P"].values
q = p * len(p) / (np.argsort(np.argsort(p)) + 1)
q = np.minimum.accumulate(q[::-1])[::-1]
res["FDR_BH"] = np.minimum(q, 1.0)
res.to_csv(os.path.join(RESD, "geneset_competitive_results.csv"), index=False)

# per-gene results for hub genes + LGALS3
hubs = ["LPCAT2", "S100A10", "CD9", "SYNGR2", "LGALS3"]
hub = gp[gp["symbol"].isin(hubs)].copy()
hub["bonf_all_genes"] = (hub["P"] * len(gp)).clip(upper=1.0)
hub.to_csv(os.path.join(RESD, "hub_genes_genelevel_results.csv"), index=False)
print(hub.to_string(index=False))

gp.sort_values("P").to_csv(os.path.join(RESD, "all_genes_results.csv.gz"), index=False, compression="gzip")
print("top genes genome-wide:")
print(gp.nsmallest(10, "P")[["symbol","chr","nSNP","P"]].to_string(index=False))
