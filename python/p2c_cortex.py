"""Phase 2c: MR + coloc for LPCAT2, CD9, SYNGR2, S100A10 using MetaBrain CORTEX
eQTL (N=2693) vs ILAE 2018 all-epilepsy GWAS (GCST007343).
Coloc restricted to +/-500 kb around gene TSS (GenePos, hg38)."""
import sys, os, gzip, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
from mr_coloc_utils import ivw, mr_egger, wald_ratio, coloc_abf, clump

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
RES = os.path.join(BASE, "results", "42_MR_coloc")

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
    g["beta_gwas"] = g["z"] * g["se_gwas"]
    return g

def harmonize(exp, gwas):
    m = exp.merge(gwas, on="SNP")
    m["EA"] = m["EA"].str.upper()
    same = m["EA"] == m["A1"]; flip = m["EA"] == m["A2"]
    m = m[same | flip].copy()
    m["beta_gwas_h"] = np.where(m["EA"] == m["A1"], m["beta_gwas"], -m["beta_gwas"])
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    other = np.where(m["EA"] == m["A1"], m["A2"], m["A1"])
    m["palindromic"] = [comp.get(e, "") == o for e, o in zip(m["EA"], other)]
    maf = np.minimum(m["eaf"], 1 - m["eaf"])
    m = m[~(m["palindromic"] & (maf > 0.42))].copy()
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
            return c, "P<5e-8 + clump r2<0.1/250kb" + (" [RELAXED: only <3 at 5e-8]" if len(c) < 3 else "")
    rel = m[m["p"] < 1e-5]
    if len(rel) == 0:
        return rel, "no SNP P<1e-5 -> insufficient instruments"
    c = clump(rel.rename(columns={"p": "p_exposure"}))
    return c, "RELAXED P<1e-5 + clump r2<0.1/250kb"

GENES = ["LPCAT2", "CD9", "SYNGR2", "S100A10"]
all_res = []
for gene in GENES:
    f = os.path.join(RES, f"exposure_{gene}_MetaBrain_cortex_raw.tsv")
    df = pd.read_csv(f, sep="\t")
    tss = int(df["GenePos"].iloc[0])
    rows = []
    for _, r in df.iterrows():
        snp = r["SNP"].split(":")
        if len(snp) >= 3 and snp[2].startswith("rs"):
            rows.append(dict(SNP=snp[2], chr=int(r["SNPChr"]), pos=int(r["SNPPos"]),
                             EA=r["SNPEffectAllele"], eaf=float(r["SNPEffectAlleleFreq"]),
                             beta=float(r["MetaBeta"]), se=float(r["MetaSE"]), p=float(r["MetaP"])))
    ex = pd.DataFrame(rows).drop_duplicates(subset=["SNP"])
    ex.to_csv(os.path.join(RES, f"exposure_{gene}_MetaBrain_cortex.csv"), index=False)
    gw = load_gwas(ex["SNP"].unique())
    m = harmonize(ex, gw)
    m.to_csv(os.path.join(RES, f"harmonized_{gene}_cortex_vs_ILAEall.csv"), index=False)
    inst, mode = pick_instruments(m)
    mr = run_mr(inst) if len(inst) else {}
    # coloc on +/-500 kb around TSS
    mc = m[(m["pos"] - tss).abs() <= 500000]
    co = coloc_abf(mc["beta"].values, mc["se"].values, mc["beta_gwas_h"].values, mc["se_gwas"].values) if len(mc) else {}
    print(f"{gene}: TSS={tss} harm={len(m)} colocSNPs(500kb)={len(mc)} instr={len(inst)} [{mode}]", flush=True)
    for k, v in mr.items():
        print(f"   {k}: b={v['b']:.4f} se={v['se']:.4f} p={v['p']:.4g}", flush=True)
    print(f"   coloc: PP.H3={co.get('PP_H3'):.4f} PP.H4={co.get('PP_H4'):.4f}", flush=True)
    all_res.append(dict(gene=gene, tissue="MetaBrain cortex EUR (N=2693)", tss_hg38=tss,
                        n_harmonized=len(m), n_coloc_snps_500kb=len(mc),
                        n_instruments=len(inst), instrument_mode=mode,
                        instruments=inst["SNP"].tolist() if len(inst) else [],
                        **{f"MR:{k}": v for k, v in mr.items()}, coloc=co))

with open(os.path.join(RES, "hub_genes_cortex_results.json"), "w") as f:
    json.dump(all_res, f, indent=1, default=float)
print("saved hub_genes_cortex_results.json")
