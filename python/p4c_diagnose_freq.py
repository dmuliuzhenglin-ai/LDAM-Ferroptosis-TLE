# -*- coding: utf-8 -*-
"""Diagnose SMR allele-frequency check failures: compare eQTL freq, GWAS freq,
and 1000G EUR genotype freq for a sample of SNPs (LGALS3 cortex)."""
import os, math, random
import pandas as pd

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
OUT = os.path.join(BASE, "results", "45_SMR_HEIDI")

eq = pd.read_csv(os.path.join(OUT, "eqtl_LGALS3_cortex.txt"), sep="\t")
ma = pd.read_csv(os.path.join(OUT, "ILAE_all_epilepsy.ma"), sep="\t")
ma = ma.set_index("SNP")

# read subset bim
bim = []
with open(os.path.join(OUT, "g1000_smr_region.bim")) as f:
    for i, line in enumerate(f):
        p = line.split("\t")
        if len(p) < 6: p = line.split()
        bim.append((i, p[1], p[4], p[5]))
nind = 503; block = (nind + 3) // 4

random.seed(1)
sample = random.sample(bim, 30)

def geno_freq(fin, idx, a1, a2):
    fin.seek(3 + idx * block)
    data = fin.read(block)
    counts = {a1: 0, a2: 0}
    n = 0
    for byte in data:
        for shift in (0, 2, 4, 6):
            g = (byte >> shift) & 3
            if g == 0: counts[a2 if False else a1] += 0; counts[a1] += 2  # hom A1
            elif g == 2: counts[a1] += 1; counts[a2] += 1  # het
            elif g == 3: counts[a2] += 2  # hom A2
            elif g == 1: n -= 1  # missing
            n += 1
            if n >= nind: break
        if n >= nind: break
    tot = counts[a1] + counts[a2]
    return counts[a1] / tot if tot else float("nan")

with open(os.path.join(OUT, "g1000_smr_region.bed"), "rb") as fin:
    rows = []
    for idx, snp, bA1, bA2 in sample:
        e = eq[eq["SNP"] == snp]
        if len(e) == 0 or snp not in ma.index: continue
        e = e.iloc[0]; g = ma.loc[snp]
        f_ref_A1 = geno_freq(fin, idx, bA1, bA2)
        rows.append((snp, bA1, bA2, round(f_ref_A1, 3),
                     e["A1"], e["A2"], round(e["Freq"], 3),
                     g["A1"], g["A2"], round(g["freq"], 3)))

print("SNP | bimA1 bimA2 f_ref(A1) | eqA1 eqA2 f_eq | gwasA1 gwasA2 f_gwas")
for r in rows:
    print(r)
