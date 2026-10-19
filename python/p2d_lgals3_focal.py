"""Phase 2d (optional): LGALS3 cortex + hippocampus eQTL vs ILAE 2018 FOCAL epilepsy (GCST007352)."""
import sys, os, gzip, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from mr_coloc_utils import ivw, mr_egger, wald_ratio, coloc_abf, clump

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
RES = os.path.join(BASE, "results", "42_MR_coloc")
TSS = 55124110  # LGALS3 TSS hg38 (MetaBrain)

def load_gwas(path, rsids):
    keep = set(rsids); out = []
    with gzip.open(path, "rt") as f:
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
    g["beta_gwas"] = g["z"] * g["se_gwas"]
    return g

def harmonize(exp, gwas):
    m = exp.merge(gwas, on="SNP")
    m["EA"] = m["EA"].str.upper()
    m = m[(m["EA"] == m["A1"]) | (m["EA"] == m["A2"])].copy()
    m["beta_gwas_h"] = np.where(m["EA"] == m["A1"], m["beta_gwas"], -m["beta_gwas"])
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    other = np.where(m["EA"] == m["A1"], m["A2"], m["A1"])
    pal = [comp.get(e, "") == o for e, o in zip(m["EA"], other)]
    maf = np.minimum(m["eaf"], 1 - m["eaf"])
    m = m[~(pd.Series(pal, index=m.index) & (maf > 0.42))].copy()
    return m

def run_mr(m):
    bx, sx = m["beta"].values, m["se"].values
    by, sy = m["beta_gwas_h"].values, m["se_gwas"].values
    out = {}
    if len(m) == 1:
        b, se, p = wald_ratio(bx[0], sx[0], by[0], sy[0])
        out["Wald ratio"] = dict(b=b, se=se, p=p)
    else:
        out["IVW (multiplicative RE)"] = ivw(bx, sx, by, sy)
        if len(m) >= 3:
            out["MR-Egger"] = mr_egger(bx, sx, by, sy)
    return out

def pick_instruments(m):
    top = m[m["p"] < 5e-8]
    if len(top) >= 1:
        c = clump(top.rename(columns={"p": "p_exposure"}))
        if len(c) >= 1:
            return c, "P<5e-8 + clump r2<0.1/250kb"
    rel = m[m["p"] < 1e-5]
    if len(rel) == 0:
        return rel, "no SNP P<1e-5 -> insufficient instruments"
    return clump(rel.rename(columns={"p": "p_exposure"})), "RELAXED P<1e-5 + clump r2<0.1/250kb"

GWAS_FOCAL = os.path.join(BASE, "data", "gwas", "ILAE2018_focal_epilepsy_METAL.gz")
res = []
for tag, expf, tissue in [
    ("cortex", "exposure_LGALS3_MetaBrain_cortex.csv", "MetaBrain cortex EUR (N=2693)"),
    ("hippocampus", "exposure_LGALS3_MetaBrain_hippocampus.csv", "MetaBrain hippocampus EUR (N=168)")]:
    ex = pd.read_csv(os.path.join(RES, expf))
    gw = load_gwas(GWAS_FOCAL, ex["SNP"].unique())
    m = harmonize(ex, gw)
    m.to_csv(os.path.join(RES, f"harmonized_LGALS3_{tag}_vs_ILAEfocal.csv"), index=False)
    inst, mode = pick_instruments(m)
    mr = run_mr(inst) if len(inst) else {}
    mc = m[(m["pos"] - TSS).abs() <= 500000]
    co = coloc_abf(mc["beta"].values, mc["se"].values, mc["beta_gwas_h"].values, mc["se_gwas"].values) if len(mc) else {}
    print(f"LGALS3 {tag} vs FOCAL: harm={len(m)} colocSNPs={len(mc)} instr={len(inst)} [{mode}]", flush=True)
    for k, v in mr.items():
        print(f"   {k}: b={v['b']:.4f} se={v['se']:.4f} p={v['p']:.4g}", flush=True)
    print(f"   coloc: PP.H3={co.get('PP_H3'):.4f} PP.H4={co.get('PP_H4'):.4f}", flush=True)
    res.append(dict(gene="LGALS3", outcome="ILAE 2018 focal epilepsy (GCST007352)", tissue=tissue,
                    n_harmonized=len(m), n_coloc_snps_500kb=len(mc), n_instruments=len(inst),
                    instrument_mode=mode, instruments=inst["SNP"].tolist() if len(inst) else [],
                    **{f"MR:{k}": v for k, v in mr.items()}, coloc=co))

with open(os.path.join(RES, "LGALS3_focal_results.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
print("saved LGALS3_focal_results.json")
