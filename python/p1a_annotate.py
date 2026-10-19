"""Phase 1 step A: build hg19 gene locations from refGene, parse ILAE GWAS,
assign SNPs to gene windows (35kb up / 10kb down, MAGMA convention)."""
import pandas as pd, numpy as np, os, sys

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
OUT = os.path.join(BASE, "data", "gwas")
os.makedirs(OUT, exist_ok=True)

# --- 1. gene locations from refGene (hg19) ---
cols = ["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd",
        "exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
rg = pd.read_csv(os.path.join(BASE, "data", "refGene_hg19.txt.gz"), sep="\t",
                 header=None, names=cols, compression="gzip")
rg["chr"] = rg["chrom"].str.replace("chr", "", regex=False)
rg = rg[rg["chr"].isin([str(i) for i in range(1, 23)])]
genes = (rg.groupby(["name2", "chr"])
           .agg(start=("txStart", "min"),
                end=("txEnd", "max"), strand=("strand", "first"))
           .reset_index().rename(columns={"name2": "symbol"}))
genes["chr"] = genes["chr"].astype(int)
genes.to_csv(os.path.join(OUT, "genes_hg19_refGene.csv"), index=False)
print("genes:", len(genes))

# --- 2. read GWAS ---
gw = pd.read_csv(os.path.join(OUT, "ILAE2018_all_epilepsy_METAL.gz"), sep="\t",
                 usecols=["CHR","BP","MarkerName","Allele1","Allele2","Freq1","Weight","Zscore","P-value"],
                 compression="gzip", dtype={"CHR": str})
gw = gw.dropna(subset=["P-value","MarkerName","CHR","BP"])
gw = gw[gw["MarkerName"].str.startswith("rs")]
gw = gw[(gw["P-value"] > 0) & (gw["P-value"] <= 1)]
gw["CHR"] = gw["CHR"].astype(int)
gw = gw[gw["CHR"].between(1, 22)]
gw = gw.drop_duplicates(subset=["MarkerName"])
print("gwas snps:", len(gw))

# --- 3. assign SNPs to gene windows ---
UP, DOWN = 35000, 10000
recs = []
for chrom, gsub in genes.groupby("chr"):
    ssub = gw[gw["CHR"] == chrom]
    if len(ssub) == 0: continue
    bp = ssub["BP"].values
    order = np.argsort(bp)
    bp_s = bp[order]
    idx = ssub.index.values[order]
    for _, g in gsub.iterrows():
        lo = np.searchsorted(bp_s, g["start"] - UP)
        hi = np.searchsorted(bp_s, g["end"] + DOWN, side="right")
        if hi > lo:
            sel = ssub.loc[idx[lo:hi], ["MarkerName","CHR","BP","Weight","Zscore","P-value","Allele1","Allele2","Freq1"]]
            sel = sel.assign(symbol=g.symbol)
            recs.append(sel)
sg = pd.concat(recs, ignore_index=True)
sg.to_csv(os.path.join(OUT, "snps_in_genes.csv"), index=False)
need = sg["MarkerName"].unique()
pd.Series(need).to_csv(os.path.join(OUT, "needed_rsids.txt"), index=False, header=False)
print("snp-gene pairs:", len(sg), "unique needed rsids:", len(need))
print("genes with >=1 SNP:", sg["symbol"].nunique())
