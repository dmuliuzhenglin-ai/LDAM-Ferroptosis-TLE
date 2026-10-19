"""Phase 1 step B: scan g1000_eur.bim for needed rsids, save row indices."""
import numpy as np, pandas as pd, os, time

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
OUT = os.path.join(BASE, "data", "gwas")
need = set(pd.read_csv(os.path.join(OUT, "needed_rsids.txt"), header=None)[0].astype(str))
print("needed:", len(need), flush=True)

idx, rsids, chrs, pos = [], [], [], []
t0 = time.time()
with open(os.path.join(BASE, "g1000_eur", "g1000_eur.bim"), "r", buffering=1024*1024) as f:
    for i, line in enumerate(f):
        parts = line.split(None, 3)
        if len(parts) >= 4 and parts[1] in need:
            p3 = parts[3].split(None, 1)[0]
            idx.append(i); rsids.append(parts[1])
            chrs.append(parts[0]); pos.append(p3)
        if i % 5000000 == 0 and i > 0:
            print(f"  scanned {i/1e6:.0f}M lines, found {len(idx)}, {time.time()-t0:.0f}s", flush=True)
df = pd.DataFrame({"bim_row": idx, "rsid": rsids, "chr": chrs, "pos": pos})
df["chr"] = pd.to_numeric(df["chr"], errors="coerce")
df = df.dropna(subset=["chr"]); df["chr"] = df["chr"].astype(int)
df = df[df["chr"].between(1, 22)]
df["pos"] = df["pos"].astype(int)
df.to_csv(os.path.join(OUT, "bim_match.csv"), index=False)
print("matched:", len(df), "in", round(time.time()-t0), "s")
