"""Phase 1 step C: MAGMA-style snp-wise-mean gene test per chromosome (resumable).
Statistic T = sum(z^2) over gene SNPs; null = weighted sum of chi2(1) with
weights = eigenvalues of LD matrix R (1000G EUR, 503 samples).
P-value: Liu-Tang-Zhang (2009) moment-matching approximation (as in SKAT),
using traces of R (eigen-free). Validated against chi2 when R=I."""
import numpy as np, pandas as pd, os, sys, time
from scipy.stats import chi2, ncx2, norm

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
OUTD = os.path.join(BASE, "data", "gwas")
RESD = os.path.join(BASE, "results", "41_MAGMA")
os.makedirs(RESD, exist_ok=True)
N = 503
ROWB = (N + 3) // 4  # 126 bytes per SNP

# --- validation of liu p-value ---
def liu_p(Q, c1, c2, c3, c4):
    if c2 <= 0:
        return 1.0
    s1 = c3 / c2**1.5
    s2 = c4 / c2**2
    if s1**2 > s2:
        a = 1.0 / (s1 - np.sqrt(max(s1**2 - s2, 1e-300)))
        d = s1 * a**3 - a**2
        l = a**2 - 2 * d
    else:
        a = 1.0 / np.sqrt(s2); d = 0.0; l = 1.0 / s2
    if l <= 0: l = 1e-6
    muQ, sdQ = c1, np.sqrt(2 * c2)
    muX, sdX = l + d, np.sqrt(2 * (l + 2 * d))
    Qn = (Q - muQ) / sdQ * sdX + muX
    if Qn < 0: return 1.0
    return float(ncx2.sf(Qn, l, d))

# validation: independent SNPs -> weights all 1 -> chi2(m)
for m in (1, 5, 50):
    Q = chi2.ppf(0.999, m)
    p = liu_p(Q, m, m, m, m)
    assert abs(p - 0.001) < 5e-5, (m, p)
print("liu validation OK (R=I matches chi2)", flush=True)

sg = pd.read_csv(os.path.join(OUTD, "snps_in_genes.csv"))
bm = pd.read_csv(os.path.join(OUTD, "bim_match.csv"))
# z from GWAS Zscore (sign irrelevant); fallback from p
sg["z"] = sg["Zscore"]
miss = sg["z"].isna()
sg.loc[miss, "z"] = -norm.isf(sg.loc[miss, "P-value"] / 2)
sg = sg.dropna(subset=["z"])

# decode table: plink bed codes 00->0(A1A1), 10->1, 11->2, 01->missing(-1)
lut = np.full((256, 4), -1, dtype=np.int8)
for b in range(256):
    for j in range(4):
        c = (b >> (2 * j)) & 3
        lut[b, j] = {0: 0, 2: 1, 3: 2, 1: -1}[c]

bed = np.memmap(os.path.join(BASE, "g1000_eur", "g1000_eur.bed"),
                dtype=np.uint8, mode="r")

todo = [int(x) for x in sys.argv[1:]] or list(range(1, 23))
for chrom in todo:
    outcsv = os.path.join(RESD, f"gene_p_chr{chrom}.csv")
    if os.path.exists(outcsv):
        print(f"chr{chrom} done, skip", flush=True); continue
    t0 = time.time()
    bm_c = bm[bm["chr"] == chrom]
    if len(bm_c) == 0: continue
    offs = 3 + bm_c["bim_row"].values.astype(np.int64) * ROWB
    raw = np.vstack([bed[o:o + ROWB] for o in offs])
    G = lut[raw].reshape(len(bm_c), ROWB * 4)[:, :N].astype(np.float32)
    del raw
    rs2loc = pd.Series(np.arange(len(bm_c)), index=bm_c["rsid"].values)
    # standardize per SNP (nan-aware), missing -> 0
    Gf = np.where(G < 0, np.nan, G)
    mu = np.nanmean(Gf, axis=1, keepdims=True)
    sd = np.nanstd(Gf, axis=1, keepdims=True)
    sd[sd == 0] = np.nan
    X = np.where(np.isnan(Gf), 0.0, (Gf - mu) / sd)
    keep = ~np.isnan(sd.ravel())
    del G, Gf
    sg_c = sg[sg["CHR"] == chrom]
    res = []
    skipped = []
    for sym, gdf in sg_c.groupby("symbol"):
        if len(gdf) > 20000:
            skipped.append(f"{sym}\tchr{chrom}\t{len(gdf)}")
            continue
        loc = rs2loc.reindex(gdf["MarkerName"].values)
        ok = loc.notna().values
        if ok.sum() == 0: continue
        mask = keep[loc.values[ok].astype(int)]
        li = loc.values[ok].astype(int)[mask]
        z = gdf["z"].values[ok][mask]
        if len(li) == 0: continue
        m = len(li)
        if m == 1:
            p = float(chi2.sf(z[0]**2, 1)); T = z[0]**2
        else:
            Xg = X[li]
            R = (Xg @ Xg.T) / (N - 1)
            T = float(np.sum(z * z))
            R2 = R @ R
            c1 = float(np.trace(R)); c2 = float(np.sum(R * R))
            c3 = float(np.sum(R * R2)); c4 = float(np.sum(R2 * R2))
            p = liu_p(T, c1, c2, c3, c4)
        res.append((sym, chrom, m, T, p))
    pd.DataFrame(res, columns=["symbol", "chr", "nSNP", "T", "P"]).to_csv(outcsv, index=False)
    if skipped:
        with open(os.path.join(RESD, "skipped_genes.txt"), "a") as fh:
            fh.write("\n".join(skipped) + "\n")
    print(f"chr{chrom}: {len(res)} genes, {len(skipped)} skipped, {time.time()-t0:.0f}s", flush=True)
