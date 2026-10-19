"""Task 8 (p20): summary-based TWAS (S-PrediXcan style) in Python.

GTEx v8 FUSION weights (Brain_Hippocampus / Brain_Cortex) x ILAE 2023 GWAS
(focal-HS / all / focal). LD from local 1000G EUR plink bed, read per-gene
window only (never whole-genome).

Statistic per gene: Z = w' z_gwas / sqrt(w' L w), L = genotype correlation
matrix of the gene's SNPs (+1e-6 ridge). BH FDR per tissue x outcome.

Memory strategy: bim and each GWAS are read once, chunked, keeping only SNPs
inside any gene window (P0-500kb .. P1+500kb); genotypes are read from the
memory-mapped bed per gene, SNP subset only.

Usage:
  gsmap_env/Scripts/python.exe scripts/p20_twas.py --tissue Hippocampus \
      [--chrom 22] [--outcomes fhs all focal]

Outputs: results/61_twas/TWAS_<tissue>_<outcome>.csv + run_<tissue>.json
"""
import os, json, time, argparse
import numpy as np
import pandas as pd
from scipy import stats as st
import rdata
from bed_reader import open_bed

BASE = os.getcwd()
OUTDIR = os.path.join(BASE, "results", "61_twas")
os.makedirs(OUTDIR, exist_ok=True)

TWAS_DIR = os.path.join(BASE, "data", "twas")
GWAS_DIR = os.path.join(BASE, "data", "gwas", "ILAE2023")
BED = os.path.join(BASE, "1000G_EUR", "g1000_eur.bed")
BIM = os.path.join(BASE, "1000G_EUR", "g1000_eur.bim")

TISSUE_DIR = {
    "Hippocampus": "GTExv8.EUR.Brain_Hippocampus",
    "Cortex": "GTExv8.EUR.Brain_Cortex",
}
OUTCOMES = {
    "fhs": "ILAE3_focal_HS_final.tbl",
    "all": "ILAE3_Caucasian_all_epilepsy_final.tbl",
    "focal": "ILAE3_Caucasian_focal_epilepsy_final.tbl",
}

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PAL = {"AT", "TA", "CG", "GC"}
WIN = 500_000

def flip(a):
    return "".join(COMP[x] for x in a)

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

# ---- GWAS parser (verified parser copied from scripts/p18_axis_cismr.py) ----
def parse_outcome_df(raw):
    raw.columns = ["CHR", "BP", "rsid", "ea", "nea", "Freq1", "FreqSE",
                   "Effective_N", "Z", "P", "tail"]
    parts = raw["tail"].str.strip().str.split(r"\s+", expand=True)
    raw["beta"] = pd.to_numeric(parts[1], errors="coerce")
    raw["se"] = pd.to_numeric(parts[2], errors="coerce")
    raw["ea"] = raw["ea"].str.upper(); raw["nea"] = raw["nea"].str.upper()
    return raw[["CHR", "BP", "rsid", "ea", "nea", "beta", "se"]]

def load_outcome_windowed(path, windows):
    """Chunked read; keep only rows whose (CHR,BP) falls in any window."""
    kept = []
    win_by_chr = {c: w[["lo", "hi"]].values for c, w in windows.groupby("chr")}
    for chunk in pd.read_csv(path, sep="\t", header=None, skiprows=1,
                             dtype=str, chunksize=1_000_000):
        df = parse_outcome_df(chunk)
        df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
        df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
        df = df.dropna(subset=["rsid", "ea", "nea", "beta", "se", "CHR", "BP"])
        mask = np.zeros(len(df), dtype=bool)
        for c, wr in win_by_chr.items():
            sub = df.CHR.values == c
            if sub.any():
                bp = df.BP.values[sub]
                inw = np.zeros(sub.sum(), dtype=bool)
                for lo, hi in wr:
                    inw |= (bp >= lo) & (bp <= hi)
                mask[sub] = inw
        kept.append(df[mask])
    out = pd.concat(kept, ignore_index=True)
    out["CHR"] = out["CHR"].astype(int)
    return out

# ---- weight model selection ----
def select_model(cvp):
    """cvp: xarray (rsq,pval) x models. Return (model, rsq, p).
    Best non-top1 model by cv R2; if its cv p>0.05 fall back to top1."""
    models = list(cvp.coords["dim_1"].values)
    rsq = cvp.sel(dim_0="rsq").values
    pv = cvp.sel(dim_0="pval").values
    non_top1 = [i for i, m in enumerate(models) if m != "top1"]
    if non_top1:
        best = non_top1[int(np.argmax(rsq[non_top1]))]
        if pv[best] <= 0.05:
            return models[best], float(rsq[best]), float(pv[best])
    i = models.index("top1")
    return "top1", float(rsq[i]), float(pv[i])

# ---- allele harmonization against weight snps (V5 = effect allele) ----
def match_sign(w5, w6, x, y):
    """Sign such that source allele x corresponds to w5 (+1) or w6 (-1);
    allows strand flip; None if no match."""
    for f in (lambda a: a, flip):
        a5, a6 = f(w5), f(w6)
        if x == a5 and y == a6:
            return 1
        if x == a6 and y == a5:
            return -1
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True, choices=list(TISSUE_DIR))
    ap.add_argument("--chrom", default=None, help="restrict to one chromosome (smoke test)")
    ap.add_argument("--outcomes", nargs="+", default=["fhs", "all", "focal"])
    args = ap.parse_args()

    tdir = TISSUE_DIR[args.tissue]
    pos = pd.read_csv(os.path.join(TWAS_DIR, f"{tdir}.pos"), sep="\t", dtype=str)
    pos["CHR"] = pos["CHR"].str.replace("chr", "", regex=False)
    pos = pos[pos.CHR.str.fullmatch(r"\d+")]
    pos["CHR"] = pos["CHR"].astype(int)
    pos["P0"] = pos["P0"].astype(int); pos["P1"] = pos["P1"].astype(int)
    if args.chrom:
        pos = pos[pos.CHR == int(args.chrom)]
    pos = pos.sort_values(["CHR", "P0"]).reset_index(drop=True)
    log(f"{args.tissue}: {len(pos)} genes" + (f" on chr{args.chrom}" if args.chrom else ""))

    # gene windows for prefiltering
    windows = pd.DataFrame({"chr": pos.CHR, "lo": pos.P0 - WIN, "hi": pos.P1 + WIN})
    win_by_chr = {c: w[["lo", "hi"]].values for c, w in windows.groupby("chr")}

    def in_windows(chr_arr, bp_arr):
        mask = np.zeros(len(chr_arr), dtype=bool)
        for c, wr in win_by_chr.items():
            sub = chr_arr == c
            if sub.any():
                bp = bp_arr[sub]
                inw = np.zeros(sub.sum(), dtype=bool)
                for lo, hi in wr:
                    inw |= (bp >= lo) & (bp <= hi)
                mask[sub] = inw
        return mask

    # bim: single chunked pass, keep window SNPs, remember global row index
    log("loading bim (window-filtered) ...")
    kept = []
    offset = 0
    for chunk in pd.read_csv(BIM, sep="\t", header=None, dtype=str,
                             names=["chr", "rsid", "cm", "pos", "a1", "a2"],
                             chunksize=2_000_000):
        chunk["chr_i"] = pd.to_numeric(chunk["chr"], errors="coerce")
        chunk["pos_i"] = pd.to_numeric(chunk["pos"], errors="coerce")
        ok = ~chunk.chr_i.isna() & ~chunk.pos_i.isna()
        chunk = chunk[ok]
        mask = in_windows(chunk.chr_i.values.astype(int), chunk.pos_i.values)
        sub = chunk[mask].copy()
        sub["brow"] = np.arange(offset, offset + len(chunk))[mask]
        kept.append(sub[["chr_i", "rsid", "pos_i", "a1", "a2", "brow"]])
        offset += len(chunk)
    bim = pd.concat(kept, ignore_index=True)
    bim = bim.drop_duplicates(subset="rsid", keep="first")
    log(f"bim kept {len(bim)} / {offset} ref SNPs")

    # outcomes: single chunked pass each, window-filtered
    gwas = {}
    for o in args.outcomes:
        gwas[o] = load_outcome_windowed(os.path.join(GWAS_DIR, OUTCOMES[o]), windows)
        gwas[o] = gwas[o].drop_duplicates(subset="rsid", keep="first")
        log(f"outcome {o}: {len(gwas[o])} window rows")

    bed = open_bed(BED, count_A1=True)  # dosage = count of bim col-5 allele
    results = {o: [] for o in args.outcomes}
    fail = {"rdat": 0, "nomatch": 0, "novar": 0, "negvar": 0}
    pal_drop_total = 0

    for chrom, gchr in pos.groupby("CHR"):
        bim_chr = bim[bim.chr_i == chrom].set_index("rsid")
        gwas_chr = {o: g[g.CHR == chrom].set_index("rsid") for o, g in gwas.items()}
        log(f"chr{chrom}: {len(gchr)} genes")

        for _, grow in gchr.iterrows():
            gid = grow["ID"]
            wgt_path = os.path.join(TWAS_DIR, grow["WGT"].replace("/", os.sep))
            try:
                rd = rdata.read_rda(wgt_path)
            except Exception:
                fail["rdat"] += 1
                continue
            model, cv_r2, cv_p = select_model(rd["cv.performance"])
            wm = rd["wgt.matrix"]
            w = wm.sel(dim_1=model).values
            snps = rd["snps"].copy()
            snps.columns = ["chr", "rsid", "cm", "pos", "a1", "a2"]
            snps["w"] = w
            snps = snps[snps.w != 0]
            if len(snps) == 0:
                fail["nomatch"] += 1
                continue
            n_wgt = len(snps)
            is_pal = (snps.a1 + snps.a2).isin(PAL)
            pal_drop_total += int(is_pal.sum())
            snps = snps[~is_pal]
            rows = []
            for _, s in snps.iterrows():
                if s.rsid not in bim_chr.index:
                    continue
                b = bim_chr.loc[s.rsid]
                sg = match_sign(s.a1, s.a2, b["a1"], b["a2"])
                if sg is None:
                    continue
                rows.append({"rsid": s.rsid, "w": s.w, "gsign": sg,
                             "brow": int(b["brow"]), "w5": s.a1, "w6": s.a2})
            if len(rows) == 0:
                fail["nomatch"] += 1
                continue
            m = pd.DataFrame(rows)
            G = bed.read(index=(slice(None), m.brow.values), dtype="float32").T
            sgn = m.gsign.values.reshape(-1, 1)
            G = np.where(sgn > 0, G, 2.0 - G)  # dosage of weight effect allele
            mu = np.nanmean(G, axis=1, keepdims=True)
            inds = np.where(np.isnan(G))
            G[inds] = np.take(mu, inds[0])
            Gc = G - mu
            sd = Gc.std(axis=1)
            keep = sd > 1e-8
            if keep.sum() == 0:
                fail["novar"] += 1
                continue
            m = m[keep].reset_index(drop=True)
            Gc = Gc[keep]; sd = sd[keep]
            L = (Gc @ Gc.T) / (Gc.shape[1] * np.outer(sd, sd))
            L[np.isnan(L)] = 0.0
            L = L + 1e-6 * np.eye(L.shape[0])
            wv = m.w.values

            for o in args.outcomes:
                go = gwas_chr[o]
                sub = m[m.rsid.isin(go.index)]
                if len(sub) == 0:
                    continue
                gs = go.loc[sub.rsid]
                zsign = np.array([np.nan if s is None else float(s) for s in
                                  (match_sign(r.w5, r.w6, e, n)
                                   for r, e, n in zip(sub.itertuples(), gs.ea, gs.nea))])
                ok = ~np.isnan(zsign)
                if ok.sum() == 0:
                    continue
                idx = sub.index.values[ok]
                z = zsign[ok] * gs.beta.values[ok] / gs.se.values[ok]
                w2 = wv[idx]
                L2 = L[np.ix_(idx, idx)]
                denom = float(w2 @ L2 @ w2)
                if denom <= 0:
                    fail["negvar"] += 1
                    continue
                ztw = float(w2 @ z / np.sqrt(denom))
                results[o].append(dict(
                    gene=gid, chr=chrom, p0=int(grow.P0), p1=int(grow.P1),
                    model=model, cv_r2=cv_r2, cv_p=cv_p,
                    nsnp=int(len(idx)), nsnp_wgt=int(n_wgt),
                    Z=ztw, P=2 * st.norm.sf(abs(ztw))))

    for o in args.outcomes:
        R = pd.DataFrame(results[o])
        if len(R) == 0:
            log(f"{args.tissue} x {o}: no results")
            continue
        R = R.sort_values("P").reset_index(drop=True)
        R["FDR"] = st.false_discovery_control(R.P.values)
        fpath = os.path.join(OUTDIR, f"TWAS_{args.tissue}_{o}.csv")
        R.to_csv(fpath, index=False)
        log(f"{args.tissue} x {o}: {len(R)} genes tested, "
            f"FDR<0.05: {(R.FDR < 0.05).sum()}, min P {R.P.min():.2e} -> {fpath}")

    log(f"fail counters: {fail}; palindromic dropped: {pal_drop_total}")
    tag = args.tissue + (f"_chr{args.chrom}" if args.chrom else "")
    with open(os.path.join(OUTDIR, f"run_{tag}.json"), "w") as f:
        json.dump({"tissue": args.tissue, "chrom": args.chrom,
                   "outcomes": args.outcomes, "fail": fail,
                   "palindromic_dropped": pal_drop_total,
                   "n_tested": {o: len(results[o]) for o in args.outcomes}}, f, indent=2)
    log("done")

if __name__ == "__main__":
    main()
