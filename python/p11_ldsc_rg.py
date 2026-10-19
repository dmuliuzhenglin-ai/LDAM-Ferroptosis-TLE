"""B5: cross-trait LD score regression (Bulik-Sullivan 2015), self-contained implementation.
Reference: gsMap packaged LDSC_resource (1000G EUR Phase3 HM3 no-HLA weights = SNP-level L2).
Munging: HM3 SNPs only, aligned to 1000G EUR A1, strand-ambiguous (A/T,C/G) dropped, MAF>0.01.
Two-step weighted regression (weights 1/L2 -> heteroscedasticity-adjusted), 200-block jackknife SE.
rg = slope12 / sqrt(slope1 * slope2) (sample size cancels).
Traits: ILAE2018_all, ILAE3_focal_HS, ILAE3_focal vs AD, STROKE, LDL, CRP.
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
RES = os.path.join(BASE, "results", "50_gsMap", "gsMap_resource")
OUT = os.path.join(BASE, "results", "57_LDSC_rg")
os.makedirs(OUT, exist_ok=True)

# ---------- reference: L2 weights + alleles ----------
wfiles = sorted(glob.glob(os.path.join(RES, "LDSC_resource", "weights_hm3_no_hla", "weights.*.l2.ldscore.gz")),
                key=lambda f: int(f.split("weights.")[1].split(".")[0]))
W = pd.concat([pd.read_csv(f, sep="\t") for f in wfiles], ignore_index=True)
print("weights SNPs:", len(W), flush=True)

bim = pd.concat([pd.read_csv(f, sep="\t", header=None,
                 names=["CHR","SNP","CM","BP","A1","A2"])
                 for f in sorted(glob.glob(os.path.join(RES, "LD_Reference_Panel", "1000G_EUR_Phase3_plink", "1000G.EUR.QC.*.bim")),
                                 key=lambda f: int(f.split("QC.")[1].split(".")[0]))], ignore_index=True)
bim["A1"] = bim["A1"].str.upper(); bim["A2"] = bim["A2"].str.upper()
ref = W.merge(bim[["SNP","A1","A2"]], on="SNP", how="inner")
ambig = ref["A1"].isin(["A","T"]) & ref["A2"].isin(["A","T"]) | ref["A1"].isin(["C","G"]) & ref["A2"].isin(["C","G"])
ref = ref[~ambig]
ref = ref[ref["MAF"] > 0.01]
print("ref after ambig/maf filter:", len(ref), flush=True)

# ---------- trait loaders: return DataFrame SNP,Z ----------
def load_metal(path, ncol, zcol):
    gw = pd.read_csv(path, sep="\t", compression="gzip" if path.endswith(".gz") else None,
                     dtype={"CHR": str})
    gw.columns = [c.strip() for c in gw.columns]
    gw = gw[["MarkerName", ncol, zcol]].rename(columns={"MarkerName": "SNP", ncol: "N", zcol: "Z"})
    return gw.dropna(subset=["SNP","Z"]).drop_duplicates("SNP")

def load_ilae3_tbl(path):
    hdr = pd.read_csv(path, sep="\t", nrows=0)
    names = [c.strip() for c in hdr.columns]
    idx = list(range(min(10, len(names))))
    gw = pd.read_csv(path, sep="\t", usecols=idx)
    gw.columns = [names[i] for i in idx]
    gw = gw[["MarkerName","Effective_N","Z-score"]].rename(
        columns={"MarkerName":"SNP","Effective_N":"N","Z-score":"Z"})
    return gw.dropna(subset=["SNP","Z"]).drop_duplicates("SNP")

def load_ebi_harm(path):
    hdr = pd.read_csv(path, sep="\t", compression="gzip", nrows=0)
    keep = [c for c in hdr.columns if c.lower() in
            ("hm_rsid","rsid","hm_beta","beta","standard_error","hm_standard_error","se",
             "p_value","hm_p_value","hm_effect_allele","effect_allele","hm_odds_ratio")]
    gw = pd.read_csv(path, sep="\t", compression="gzip", usecols=keep, low_memory=False)
    cols = {c.lower(): c for c in gw.columns}
    rsid = cols.get("hm_rsid") or cols.get("rsid")
    beta = cols.get("hm_beta") or cols.get("beta")
    se = cols.get("standard_error") or cols.get("hm_standard_error") or cols.get("se")
    p = cols.get("p_value") or cols.get("hm_p_value")
    ea = cols.get("hm_effect_allele") or cols.get("effect_allele")
    out = pd.DataFrame({"SNP": gw[rsid]})
    if beta and se:
        b = pd.to_numeric(gw[beta], errors="coerce")
        s = pd.to_numeric(gw[se], errors="coerce")
        out["Z"] = b / s
    elif p:
        pv = pd.to_numeric(gw[p], errors="coerce")
        out["Z"] = -stats.norm.isf(pv / 2) * np.sign(pd.to_numeric(gw.get(cols.get("hm_odds_ratio")), errors="coerce").apply(np.log) if cols.get("hm_odds_ratio") else 1)
    out["EA"] = gw[ea].str.upper() if ea else None
    return out.dropna(subset=["SNP","Z"]).drop_duplicates("SNP")

def load_with_alleles_metal(path):
    gw = pd.read_csv(path, sep="\t", compression="gzip" if path.endswith(".gz") else None, dtype={"CHR": str})
    gw.columns = [c.strip() for c in gw.columns]
    ncol = "Weight" if "Weight" in gw.columns else "Effective_N"
    zcol = "Zscore" if "Zscore" in gw.columns else "Z-score"
    out = pd.DataFrame({"SNP": gw["MarkerName"], "Z": gw[zcol], "EA": gw["Allele1"].str.upper(),
                        "OA": gw["Allele2"].str.upper(), "N": gw[ncol]})
    return out.dropna(subset=["SNP","Z"]).drop_duplicates("SNP")

D = os.path.join(BASE, "data")
TRAITS = {
    "ILAE2018_all":  lambda: load_with_alleles_metal(os.path.join(D, "gwas", "ILAE2018_all_epilepsy_METAL.gz")),
    "ILAE3_all":     lambda: load_with_alleles_metal(os.path.join(D, "gwas", "ILAE2023", "ILAE3_Caucasian_all_epilepsy_final.tbl")),
    "ILAE3_focal":   lambda: load_with_alleles_metal(os.path.join(D, "gwas", "ILAE2023", "ILAE3_Caucasian_focal_epilepsy_final.tbl")),
    "ILAE3_focalHS": lambda: load_with_alleles_metal(os.path.join(D, "gwas", "ILAE2023", "ILAE3_focal_HS_final.tbl")),
    "AD":     lambda: load_ebi_harm(os.path.join(D, "ldsc", "AD_GCST90027158.h.tsv.gz")),
    "STROKE": lambda: load_ebi_harm(os.path.join(D, "ldsc", "STROKE_GCST005838.h.tsv.gz")),
    "LDL":    lambda: load_ebi_harm(os.path.join(D, "ldsc", "LDL_GCST90092883.h.tsv.gz")),
    "CRP":    lambda: load_ebi_harm(os.path.join(D, "ldsc", "CRP_GCST90029070.h.tsv.gz")),
}

# ---------- munge onto reference ----------
def munge(df):
    df = df.copy()
    df["EA"] = df["EA"].astype(str).str.upper()
    m = ref.merge(df[["SNP","Z","EA"]], on="SNP", how="inner")
    flip = m["EA"] == m["A2"]
    keep = (m["EA"] == m["A1"]) | flip
    m = m[keep]
    m.loc[flip[keep].values, "Z"] = -m.loc[flip[keep].values, "Z"]
    return m[["SNP","L2","MAF","Z"]].reset_index(drop=True)

# ---------- LDSC core ----------
def wls(x, y, w):
    sw = w.sum()
    xm = (w * x).sum() / sw; ym = (w * y).sum() / sw
    b = (w * (x - xm) * (y - ym)).sum() / (w * (x - xm) ** 2).sum()
    a = ym - b * xm
    return a, b

def jackknife_fixed(x, y, w, nblocks=200):
    """WLS with fixed weights; delete-block jackknife SE of slope."""
    a0, b0 = wls(x, y, w)
    n = len(x)
    blocks = np.array_split(np.arange(n), nblocks)
    ests = []
    for bl in blocks:
        mask = np.ones(n, bool); mask[bl] = False
        ests.append(wls(x[mask], y[mask], w[mask])[1])
    ests = np.array(ests)
    se = np.sqrt((nblocks - 1) / nblocks * np.sum((ests - ests.mean()) ** 2))
    return a0, b0, se, ests

def h2_est(z, l2, w_l2):
    x, y = l2, z ** 2
    a, b = wls(x, y, 1.0 / w_l2)                      # step 1
    het = np.maximum(2 * (1 + b * x) ** 2, 1e-3)      # approx var(chi2 | L)
    w = 1.0 / (w_l2 * het)
    a, b, se, _ = jackknife_fixed(x, y, w)
    return a, b, se

def rg_est_jk(z1, z2, l2, b1, b2):
    x, y = l2, z1 * z2
    het = np.maximum((1 + b1 * x) * (1 + b2 * x), 1e-3)
    w = 1.0 / (l2 * het)
    a, b, se, ests = jackknife_fixed(x, y, w)
    return a, b, se, ests

# ---------- run ----------
cache = {}
def get_trait(name):
    if name not in cache:
        f = os.path.join(OUT, f"munged_{name}.csv.gz")
        if os.path.exists(f):
            cache[name] = pd.read_csv(f)
        else:
            df = TRAITS[name]()
            cache[name] = munge(df)
            cache[name].to_csv(f, index=False, compression="gzip")
        print(name, "munged SNPs:", len(cache[name]), flush=True)
    return cache[name]

FOCAL = ["ILAE2018_all", "ILAE3_all", "ILAE3_focal", "ILAE3_focalHS"]
EXT = ["AD", "STROKE", "LDL", "CRP"]

h2tab = {}
for t in FOCAL + EXT:
    try:
        d = get_trait(t)
    except Exception as e:
        print(t, "LOAD FAIL:", e, flush=True); continue
    a, b, se = h2_est(d["Z"].values, d["L2"].values, d["L2"].values)
    zsc = b / se if se > 0 else np.nan
    h2tab[t] = dict(trait=t, nsnp=len(d), slope=b, slope_se=se, slope_z=zsc, intercept=a)
    print(t, h2tab[t], flush=True)
pd.DataFrame(h2tab.values()).to_csv(os.path.join(OUT, "h2_slopes.csv"), index=False)

rows = []
import itertools
PAIRS = [(a, b) for a, b in itertools.combinations(FOCAL, 2)] + \
        [(tf, te) for tf in FOCAL for te in EXT]
for tf, te in PAIRS:
        if tf not in h2tab or te not in h2tab:
            continue
        d1, d2 = get_trait(tf), get_trait(te)
        m = d1.merge(d2, on="SNP", suffixes=("_1", "_2"))
        x = m["L2_1"].values
        b1f, b2f = h2tab[tf]["slope"], h2tab[te]["slope"]
        a, b, se_b, ests = rg_est_jk(m["Z_1"].values, m["Z_2"].values, x, b1f, b2f)
        denom = np.sqrt(b1f * b2f)
        rg = b / denom
        rg_ests = ests / denom
        se_rg = np.sqrt((199 / 200) * np.sum((rg_ests - rg_ests.mean()) ** 2))
        p = 2 * stats.norm.sf(abs(rg / se_rg)) if se_rg > 0 else np.nan
        rows.append(dict(trait1=tf, trait2=te, nsnp=len(m), rg=rg, se=se_rg, z=rg / se_rg, p=p,
                         gcov_slope=b, gcov_slope_se=se_b))
        print(rows[-1], flush=True)
res = pd.DataFrame(rows)
if len(res):
    from statsmodels.stats.multitest import multipletests
    res["fdr"] = multipletests(res["p"], method="fdr_bh")[1]
res.to_csv(os.path.join(OUT, "ldsc_rg_results.csv"), index=False)
print("B5_DONE", flush=True)
