# -*- coding: utf-8 -*-
"""Subset g1000_eur bed/bim/fam to the SNPs needed for SMR HEIDI (region SNPs),
so the SMR Windows binary does not have to scan a 658MB bim / 2.8GB bed."""
import os, math

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
REF = os.path.join(BASE, "g1000_eur", "g1000_eur")
OUT = os.path.join(BASE, "results", "45_SMR_HEIDI")
MA = os.path.join(OUT, "ILAE_all_epilepsy.ma")

# SNPs needed = union of exposure SNPs (same set used for .ma)
need = set()
with open(MA) as f:
    next(f)
    for line in f:
        need.add(line.split("\t", 1)[0])
print("needed:", len(need))

fam = open(REF + ".fam").read().strip().split("\n")
nind = len(fam)
block = (nind + 3) // 4
print("individuals:", nind, "block bytes:", block)

keep_idx = []   # (idx, bim_line)
with open(REF + ".bim") as f:
    for i, line in enumerate(f):
        parts = line.split("\t")
        if len(parts) < 6:
            parts = line.split()
        if parts[1] in need:
            keep_idx.append((i, line if line.endswith("\n") else line + "\n"))
print("matched in bim:", len(keep_idx))

with open(REF + ".bed", "rb") as fin, \
     open(os.path.join(OUT, "g1000_smr_region.bed"), "wb") as fout, \
     open(os.path.join(OUT, "g1000_smr_region.bim"), "w") as fbim:
    magic = fin.read(3)
    fout.write(magic)
    for idx, line in keep_idx:
        fin.seek(3 + idx * block)
        fout.write(fin.read(block))
        fbim.write(line)

import shutil
shutil.copy(REF + ".fam", os.path.join(OUT, "g1000_smr_region.fam"))
print("done ->", os.path.join(OUT, "g1000_smr_region.{bed,bim,fam}"))
