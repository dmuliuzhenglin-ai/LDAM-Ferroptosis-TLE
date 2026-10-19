"""Phase 2: LGALS3 MR + coloc, exposure=MetaBrain hippocampus (local chr14) and
cortex (local SMR_analysis txt), outcome=ILAE 2018 all-epilepsy GWAS (local)."""
import sys, os, gzip, io
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
from mr_coloc_utils import ivw, mr_egger, wald_ratio, coloc_abf, clump, _coloc_selftest

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
RES = os.path.join(BASE, "results", "42_MR_coloc")
os.makedirs(RES, exist_ok=True)
assert _coloc_selftest()
print("coloc self-test OK", flush=True)

# ---------- load GWAS (only needed rsids) ----------
def load_gwas(rsids):
    keep = set(rsids)
    out = []
    with gzip.open(os.path.join(BASE, "data", "gwas", "ILAE2018_all_epilepsy_METAL.gz"), "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        iM = header.index("MarkerName"); iA1 = header.index("Allele1"); iA2 = header.index("Allele2")
        iF = header.index("Freq1"); iW = header.index("Weight"); iZ = header.index("Zscore")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[iM] in keep:
                out.append((p[iM], p[iA1].upper(), p[iA2].upper(),
                            float(p[iF]), float(p[iW]), float(p[iZ])))
    g = pd.DataFrame(out, columns=["SNP", "A1", "A2", "freq1", "N", "z"])
    g["se_gwas"] = 1.0 / np.sqrt(2 * g["freq1"] * (1 - g["freq1"]) * (g["N"] + g["z"] ** 2))
    g["beta_gwas"] = g["z"] * g["se_gwas"]  # effect w.r.t. A1
    return g

# ---------- harmonize ----------
def harmonize(exp, gwas, ea_col, eaf_col):
    m = exp.merge(gwas, on="SNP")
    m[ea_col] = m[ea_col].str.upper()
    m["A1"] = m["A1"].str.upper(); m["A2"] = m["A2"].str.upper()
    same = m[ea_col] == m["A1"]
    flip = m[ea_col] == m["A2"]
    m = m[same | flip].copy()
    m["beta_gwas_h"] = np.where(m[ea_col] == m["A1"], m["beta_gwas"], -m["beta_gwas"])
    pal = m[ea_col].isin(["A", "T"]) & m["A2" if False else ea_col].notna()
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    other = np.where(m[ea_col] == m["A1"], m["A2"], m["A1"])
    m["palindromic"] = [comp.get(e, "") == o for e, o in zip(m[ea_col], other)]
    maf = np.minimum(m[eaf_col], 1 - m[eaf_col])
    drop_pal = m["palindromic"] & (maf > 0.42)
    m = m[~drop_pal].copy()
    return m

def run_mr(m, bxcol, sxcol, bycol="beta_gwas_h", sycol="se_gwas"):
    bx, sx = m[bxcol].values, m[sxcol].values
    by, sy = m[bycol].values, m[sycol].values
    out = {}
    if len(m) == 1:
        b, se, p = wald_ratio(bx[0], sx[0], by[0], sy[0])
        out["Wald ratio"] = dict(b=b, se=se, p=p)
    else:
        out["IVW (multiplicative RE)"] = ivw(bx, sx, by, sy)
        if len(m) >= 3:
            out["MR-Egger"] = mr_egger(bx, sx, by, sy)
    return out

def pick_instruments(m, pcol, pos_col, chr_col):
    top = m[m[pcol] < 5e-8]
    mode = "P<5e-8, no clumping needed" if len(top) >= 3 else ""
    if len(top) >= 3:
        # still clump to ensure independence
        c = clump(top.rename(columns={pcol: "p_exposure", pos_col: "pos", chr_col: "chr"}))
        if len(c) >= 1:
            return c, "P<5e-8 + clump r2<0.1/250kb"
    rel = m[m[pcol] < 1e-5]
    if len(rel) == 0:
        return rel, "no SNP P<1e-5 -> insufficient instruments"
    c = clump(rel.rename(columns={pcol: "p_exposure", pos_col: "pos", chr_col: "chr"}))
    return c, "RELAXED P<1e-5 + clump r2<0.1/250kb"

results = []

# ========== exposure 1: MetaBrain hippocampus (local chr14, N=168) ==========
rows = []
with gzip.open(os.path.join(BASE, "2021-07-23-hippocampus-EUR-30PCs-chr14.txt.gz"), "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(header)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[idx["GeneSymbol"]] == "LGALS3":
            snp = p[idx["SNP"]].split(":")
            rows.append(dict(SNP=snp[2], chr=int(p[idx["SNPChr"]]), pos=int(p[idx["SNPPos"]]),
                             EA=p[idx["SNPEffectAllele"]], eaf=float(p[idx["SNPEffectAlleleFreq"]]),
                             beta=float(p[idx["MetaBeta"]]), se=float(p[idx["MetaSE"]]),
                             p=float(p[idx["MetaP"]])))
hs = pd.DataFrame(rows).drop_duplicates(subset=["SNP"])
hs = hs[hs["SNP"].str.startswith("rs")]
hs.to_csv(os.path.join(RES, "exposure_LGALS3_MetaBrain_hippocampus.csv"), index=False)
print("hippocampus LGALS3 SNPs:", len(hs), "minP:", hs["p"].min())

gw = load_gwas(hs["SNP"].unique())
print("GWAS matched:", len(gw))
m = harmonize(hs, gw, "EA", "eaf")
m.to_csv(os.path.join(RES, "harmonized_LGALS3_hippocampus_vs_ILAEall.csv"), index=False)
inst, mode = pick_instruments(m, "p", "pos", "chr")
print("hippocampus instruments:", len(inst), mode)
mr = run_mr(inst, "beta", "se") if len(inst) else {}
co = coloc_abf(m["beta"].values, m["se"].values, m["beta_gwas_h"].values, m["se_gwas"].values)
print("hippocampus coloc:", co)
results.append(dict(gene="LGALS3", tissue="MetaBrain hippocampus EUR (N=168)",
                    n_harmonized=len(m), n_instruments=len(inst), instrument_mode=mode,
                    **{f"MR:{k}": v for k, v in mr.items()}, coloc=co))

# ========== exposure 2: MetaBrain cortex (local SMR txt, N=2693) ==========
cx = pd.read_csv(os.path.join(BASE, "SMR_analysis", "LGALS3_cortex_eqtl.txt"), sep=r"\s+")
cx = cx.rename(columns={"SNP": "SNP", "Chr": "chr", "BP": "pos", "A1": "EA", "Freq": "eaf",
                        "b": "beta", "se": "se", "p": "p"})
cx = cx[["SNP", "chr", "pos", "EA", "eaf", "beta", "se", "p"]].drop_duplicates(subset=["SNP"])
cx.to_csv(os.path.join(RES, "exposure_LGALS3_MetaBrain_cortex.csv"), index=False)
gw2 = load_gwas(cx["SNP"].unique())
m2 = harmonize(cx, gw2, "EA", "eaf")
m2.to_csv(os.path.join(RES, "harmonized_LGALS3_cortex_vs_ILAEall.csv"), index=False)
inst2, mode2 = pick_instruments(m2, "p", "pos", "chr")
print("cortex instruments:", len(inst2), mode2)
mr2 = run_mr(inst2, "beta", "se") if len(inst2) else {}
co2 = coloc_abf(m2["beta"].values, m2["se"].values, m2["beta_gwas_h"].values, m2["se_gwas"].values)
print("cortex coloc:", co2)
results.append(dict(gene="LGALS3", tissue="MetaBrain cortex EUR (N=2693)",
                    n_harmonized=len(m2), n_instruments=len(inst2), instrument_mode=mode2,
                    **{f"MR:{k}": v for k, v in mr2.items()}, coloc=co2))

import json
with open(os.path.join(RES, "LGALS3_results.json"), "w") as f:
    json.dump(results, f, indent=1, default=float)
print(json.dumps(results, indent=1, default=float)[:3000])
