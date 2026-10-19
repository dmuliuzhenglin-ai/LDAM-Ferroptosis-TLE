"""p21: TWAS implementation validation against ILAE2023 all-epilepsy GWAS.

Logic: summary-based TWAS with correct LD + harmonization should rank genes
near known genome-wide significant GWAS loci higher than random. We:
  1. derive independent lead loci (P<5e-8, +-1Mb clump) from ILAE3 all-epilepsy
  2. for each locus, find tested TWAS genes overlapping +-500kb and record the
     best (min P) gene rank percentile genome-wide
  3. specifically check the 2q24.3 SCN1A/SCN2A/SCN3A sodium-channel cluster
  4. locate all FDR<0.05 TWAS genes relative to GWAS loci
Output: results/61_twas/validation_notes.txt (+ locus table CSV)
"""
import os
import numpy as np
import pandas as pd

BASE = os.getcwd()
OUT = os.path.join(BASE, "results", "61_twas")
GWAS = os.path.join(BASE, "data", "gwas", "ILAE2023",
                    "ILAE3_Caucasian_all_epilepsy_final.tbl")
GTF = os.path.join(BASE, "results", "50_gsMap", "gsMap_resource",
                   "genome_annotation", "gtf",
                   "gencode.v46lift37.basic.annotation.gtf")

def ens2sym():
    m = {}
    with open(GTF, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            attrs = c[8]
            gid = attrs.split('gene_id "')[1].split('"')[0].split(".")[0].split("_")[0]
            gn = attrs.split('gene_name "')[1].split('"')[0] if 'gene_name "' in attrs else gid
            m[gid] = gn
    return m

# p18 verified parser
def load_outcome(path):
    raw = pd.read_csv(path, sep="\t", header=None, skiprows=1, dtype=str)
    raw.columns = ["CHR", "BP", "rsid", "ea", "nea", "Freq1", "FreqSE",
                   "Effective_N", "Z", "P", "tail"]
    parts = raw["tail"].str.strip().str.split(r"\s+", expand=True)
    raw["beta"] = pd.to_numeric(parts[1], errors="coerce")
    raw["se"] = pd.to_numeric(parts[2], errors="coerce")
    raw["Pgwas"] = pd.to_numeric(raw["P"], errors="coerce")
    raw["CHR"] = pd.to_numeric(raw["CHR"], errors="coerce")
    raw["BP"] = pd.to_numeric(raw["BP"], errors="coerce")
    return raw.dropna(subset=["CHR", "BP", "Pgwas"])

def clump_loci(g, pthr=5e-8, win=1_000_000):
    hits = g[g.Pgwas < pthr].sort_values(["CHR", "BP"])
    loci = []
    for (c, grp) in hits.groupby("CHR"):
        cur = None
        for _, r in grp.iterrows():
            if cur is None or r.BP - cur["end"] > win:
                if cur is not None:
                    loci.append(cur)
                cur = {"chr": int(c), "start": int(r.BP), "end": int(r.BP),
                       "lead_p": r.Pgwas, "lead_pos": int(r.BP), "lead_rsid": r.rsid}
            else:
                cur["end"] = int(r.BP)
                if r.Pgwas < cur["lead_p"]:
                    cur.update(lead_p=r.Pgwas, lead_pos=int(r.BP), lead_rsid=r.rsid)
        if cur is not None:
            loci.append(cur)
    return pd.DataFrame(loci)

def validate(tissue, R, loci, s2m, lines):
    n = len(R)
    R = R.copy()
    R["rank_pct"] = R.P.rank(pct=True)  # fraction of genes with P <= this P
    lines.append(f"=== {tissue} x all-epilepsy: {n} genes tested ===")
    # per-locus best gene
    rows = []
    for _, L in loci.iterrows():
        sub = R[(R.chr == L.chr) & (R.p1 >= L.start - 500_000) & (R.p0 <= L.end + 500_000)]
        if len(sub) == 0:
            rows.append({**L.to_dict(), "n_twas_genes": 0})
            continue
        b = sub.loc[sub.P.idxmin()]
        rows.append({**L.to_dict(), "n_twas_genes": len(sub),
                     "best_gene": b.gene, "best_sym": s2m.get(b.gene.split(".")[0], b.gene),
                     "best_P": b.P, "best_Z": b.Z, "best_FDR": b.FDR,
                     "best_rank_pct": b.rank_pct})
    LT = pd.DataFrame(rows)
    tested = LT[LT.n_twas_genes > 0]
    lines.append(f"GWAS lead loci (P<5e-8, +-1Mb clump): {len(loci)}; "
                 f"with >=1 tested TWAS gene nearby: {len(tested)}")
    if len(tested):
        med = tested.best_rank_pct.median()
        frac_top10 = (tested.best_rank_pct <= 0.10).mean()
        frac_top25 = (tested.best_rank_pct <= 0.25).mean()
        lines.append(f"best-gene rank percentile at GWAS loci: median {med:.3f}; "
                     f"in top 10% genome-wide: {frac_top10:.2f} ({int(frac_top10*len(tested))}/{len(tested)}); "
                     f"in top 25%: {frac_top25:.2f} ({int(frac_top25*len(tested))}/{len(tested)})")
        lines.append("(null expectation: median 0.50, top-10% fraction 0.10, top-25% fraction 0.25)")
    # SCN1A cluster (2q24.3): SCN2A 165.2Mb, SCN3A 165.9Mb, SCN1A 166.85Mb
    cl = R[(R.chr == 2) & (R.p0 >= 164_000_000) & (R.p1 <= 169_000_000)]
    lines.append(f"\n2q24.3 SCN cluster (chr2:164-169Mb): {len(cl)} tested genes")
    for _, r in cl.sort_values("P").iterrows():
        sym = s2m.get(r.gene.split(".")[0], r.gene)
        mark = " <-- SCN gene" if sym in ("SCN1A", "SCN2A", "SCN3A", "SCN7A", "SCN9A") else ""
        lines.append(f"  {sym} ({r.gene}): Z={r.Z:+.2f}, P={r.P:.2e}, "
                     f"rank_pct={r.rank_pct:.3f}{mark}")
    # FDR-significant genes vs loci
    sig = R[R.FDR < 0.05]
    lines.append(f"\nFDR<0.05 TWAS genes: {len(sig)}")
    for _, r in sig.sort_values("P").iterrows():
        sym = s2m.get(r.gene.split(".")[0], r.gene)
        near = loci[(loci.chr == r.chr) & (loci.end + 500_000 >= r.p0) & (loci.start - 500_000 <= r.p1)]
        tag = (f"overlaps GWAS locus {near.iloc[0].chr}:{near.iloc[0].start}-{near.iloc[0].end} "
               f"(lead {near.iloc[0].lead_rsid}, Pgwas={near.iloc[0].lead_p:.1e})") if len(near) else "no GWAS locus within 500kb"
        lines.append(f"  {sym} ({r.gene}) chr{r.chr}:{r.p0}-{r.p1}: Z={r.Z:+.2f}, "
                     f"P={r.P:.2e}, FDR={r.FDR:.3f} -- {tag}")
    return LT

def main():
    s2m = ens2sym()
    print("gencode map:", len(s2m))
    g = load_outcome(GWAS)
    print("GWAS rows:", len(g))
    loci = clump_loci(g)
    print("loci:", len(loci))
    lines = ["TWAS implementation validation (p21)",
             "GWAS: ILAE 2023 all-epilepsy (EUR). Tissues: GTEx v8 hippocampus, cortex.",
             "Z_twas = w'z / sqrt(w'Lw); L = 1000G EUR per-gene correlation (+1e-6 ridge).",
             ""]
    tabs = {}
    for tissue in ["Hippocampus", "Cortex"]:
        f = os.path.join(OUT, f"TWAS_{tissue}_all.csv")
        if not os.path.exists(f):
            lines.append(f"=== {tissue}: result file missing ===")
            continue
        R = pd.read_csv(f)
        tabs[tissue] = validate(tissue, R, loci, s2m, lines)
        lines.append("")
    for tissue, LT in tabs.items():
        LT.to_csv(os.path.join(OUT, f"validation_loci_{tissue}.csv"), index=False)
    txt = "\n".join(lines)
    with open(os.path.join(OUT, "validation_notes.txt"), "w", encoding="utf-8") as f:
        f.write(txt)
    print(txt)

if __name__ == "__main__":
    main()
