# -*- coding: utf-8 -*-
"""
Task C: Spatially variable gene (SVG) statistics for axis genes.

For each sample (16 rat GSE301260 + 2 human GSE288888), compute Moran's I
(KNN k=6, symmetrised row-standardised weights, 499 permutations, two-sided)
per axis gene on log-normalised expression (same QC as p31).

Axis gene sets:
  LDAM signature   genesets/LDAM_geneset.txt
  Ferroptosis core curated core machinery list (intersected with
                   genesets/Ferroptosis_geneset.txt)
  Module5          results/03_WGCNA_LDAM_hub_module_genes.txt
  Hub              LGALS3, SPP1, APOE + results/04_hub_genes.txt

Rat symbols resolved human->rat via results/66_spatial_rat/ortholog_map.json
cache (exact/title-case/case-insensitive first, mygene.info fallback, cache
updated on disk).

Outputs: results/71_domain_dynamics_svg/parts/svg_<sample>.csv
Run: gsmap_env\\Scripts\\python.exe scripts\\p36_svg.py --species rat --start 0 --end 4
"""
import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sparse
import scanpy as sc

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results" / "71_domain_dynamics_svg" / "parts"
OUT.mkdir(parents=True, exist_ok=True)
ORTHO_CACHE = BASE / "results" / "66_spatial_rat" / "ortholog_map.json"

RAT_DIR = BASE / "data" / "spatial" / "GSE301260"
HUM_DIR = BASE / "data" / "spatial" / "GSE288888"

RAT_SAMPLES = ["GSM9078930_A_L1_S1", "GSM9078931_A_L2_S5",
               "GSM9078932_A_L3_S9", "GSM9078933_A_L4_S13",
               "GSM9078934_B_L1_S2", "GSM9078935_B_L2_S6",
               "GSM9078936_B_L3_S10", "GSM9078937_B_L4_S14",
               "GSM9078938_C_L1_S3", "GSM9078939_C_L2_S7",
               "GSM9078940_C_L3_S11", "GSM9078941_C_L4_S15",
               "GSM9078942_D_L1_S4", "GSM9078943_D2_L2_S8",
               "GSM9078944_D2_L3_S12", "GSM9078945_D_L4_S16"]
RAT_META = {  # sample -> (treatment, timepoint)
    "GSM9078930_A_L1_S1": ("SE", "D5"), "GSM9078931_A_L2_S5": ("SE", "D10"),
    "GSM9078932_A_L3_S9": ("SE", "D5"), "GSM9078933_A_L4_S13": ("SE", "D10"),
    "GSM9078934_B_L1_S2": ("CTRL", "D5"), "GSM9078935_B_L2_S6": ("CTRL", "D10"),
    "GSM9078936_B_L3_S10": ("CTRL", "D5"), "GSM9078937_B_L4_S14": ("CTRL", "D10"),
    "GSM9078938_C_L1_S3": ("SE", "D20"), "GSM9078939_C_L2_S7": ("SE", "D40"),
    "GSM9078940_C_L3_S11": ("SE", "D40"), "GSM9078941_C_L4_S15": ("SE", "D20"),
    "GSM9078942_D_L1_S4": ("CTRL", "D20"), "GSM9078943_D2_L2_S8": ("CTRL", "D40"),
    "GSM9078944_D2_L3_S12": ("CTRL", "D40"), "GSM9078945_D_L4_S16": ("CTRL", "D20"),
}
HUM_SAMPLES = {"human_CTRL": ("GSM8778056_A18-43", "CTRL"),
               "human_HS": ("GSM8778057_276A", "HS")}

FERRO_CORE = ["GPX4", "ACSL4", "TFRC", "FTH1", "FTL", "NCOA4", "HMOX1",
              "SLC7A11", "SLC3A2", "ALOX15", "ALOX12", "LPCAT3", "AIFM2",
              "DHODH", "GCH1", "NFE2L2", "KEAP1", "ATG5", "ATG7", "BECN1",
              "MAP1LC3B", "SAT1", "GLS2", "IREB2", "STEAP3", "VDAC2",
              "CYBB", "NOX4", "CHAC1", "HSPB1", "CARS1", "ACSF2", "DPP4",
              "HSPA5", "POR", "CYB5R1"]


def read_lines(p):
    return [l.strip() for l in open(p, encoding="utf-8") if l.strip()]


def build_gene_sets():
    ld = read_lines(BASE / "genesets" / "LDAM_geneset.txt")
    ferro_all = set(read_lines(BASE / "genesets" / "Ferroptosis_geneset.txt"))
    fc = [g for g in FERRO_CORE if g in ferro_all]
    fc_missing = [g for g in FERRO_CORE if g not in ferro_all]
    m5 = read_lines(BASE / "results" / "03_WGCNA_LDAM_hub_module_genes.txt")
    hub = ["LGALS3", "SPP1", "APOE"] + read_lines(BASE / "results" / "04_hub_genes.txt")
    sets = {"LDAM": ld, "Ferroptosis_core": fc, "Module5": m5, "Hub": hub}
    union = sorted(set(ld) | set(fc) | set(m5) | set(hub))
    print(f"[genes] LDAM={len(ld)} FerroCore={len(fc)} M5={len(m5)} "
          f"Hub={len(hub)} union={len(union)} ferro_core_not_in_geneset={fc_missing}")
    return sets, union


# ------------------------------------------------------ ortholog mapping
def map_orthologs(human_genes, rat_symbols):
    rat_set = set(rat_symbols)
    rat_lower = {}
    for s in rat_set:
        rat_lower.setdefault(s.lower(), s)
    cache = json.loads(ORTHO_CACHE.read_text(encoding="utf-8")) \
        if ORTHO_CACHE.exists() else {}
    mapping, unmatched = {}, []
    for g in human_genes:
        if g in cache and cache[g] and cache[g] in rat_set:
            mapping[g] = cache[g]
        elif g in rat_set:
            mapping[g] = g
        elif g.title() in rat_set:
            mapping[g] = g.title()
        elif g.lower() in rat_lower:
            mapping[g] = rat_lower[g.lower()]
        elif g in cache and cache[g] is None:
            mapping[g] = None
        else:
            unmatched.append(g)
    if unmatched:
        try:
            import requests
            found = {}
            for i in range(0, len(unmatched), 400):
                chunk = unmatched[i:i + 400]
                r = requests.post(
                    "https://mygene.info/v3/query",
                    data={"q": ",".join(chunk), "scopes": "symbol,alias",
                          "fields": "symbol,pantherdb.ortholog",
                          "species": "human"}, timeout=30)
                r.raise_for_status()
                rat_entrez = {}
                for hit in r.json():
                    q = hit.get("query")
                    orth = hit.get("pantherdb", {}).get("ortholog", [])
                    if isinstance(orth, dict):
                        orth = [orth]
                    for o in orth:
                        if o.get("species") == "RAT":
                            rat_entrez[q] = str(o.get("entrezgene"))
                            break
                if rat_entrez:
                    rr = requests.post(
                        "https://mygene.info/v3/query",
                        data={"q": ",".join(rat_entrez.values()),
                              "scopes": "entrezgene", "fields": "symbol",
                              "species": "rat"}, timeout=30)
                    rr.raise_for_status()
                    ez2sym = {h.get("query"): h.get("symbol")
                              for h in rr.json() if "symbol" in h}
                    for g, ez in rat_entrez.items():
                        sym = ez2sym.get(ez)
                        if sym and sym in rat_set:
                            found[g] = sym
                        elif sym and sym.lower() in rat_lower:
                            found[g] = rat_lower[sym.lower()]
                time.sleep(0.34)
            for g in unmatched:
                mapping[g] = found.get(g)
                cache[g] = found.get(g)
        except Exception as e:
            print(f"[warn] mygene failed: {e}")
            for g in unmatched:
                mapping[g] = None
                cache[g] = None
        ORTHO_CACHE.write_text(json.dumps(cache, indent=1), encoding="utf-8")
    n_ok = sum(1 for v in mapping.values() if v)
    print(f"[map] human->rat {n_ok}/{len(human_genes)} mapped")
    return mapping


# ------------------------------------------------------------- loading
def read_positions(path):
    with gzip.open(path, "rt") as f:
        first = f.readline()
    cols = ["barcode", "in_tissue", "array_row", "array_col",
            "pxl_row_in_fullres", "pxl_col_in_fullres"]
    if first.lower().startswith("barcode"):
        pos = pd.read_csv(path)
        pos.columns = [c.strip() for c in pos.columns]
    else:
        pos = pd.read_csv(path, header=None, names=cols)
    return pos


def load_sample(prefix, pos_file):
    with gzip.open(f"{prefix}_matrix.mtx.gz", "rb") as f:
        mat = scipy.io.mmread(f)
    X = sparse.csr_matrix(mat).T.tocsr()
    feats = pd.read_csv(f"{prefix}_features.tsv.gz", sep="\t", header=None,
                        names=["id", "symbol", "ftype"])
    bc = pd.read_csv(f"{prefix}_barcodes.tsv.gz", header=None,
                     names=["barcode"])
    pos = read_positions(f"{prefix}_{pos_file}")
    ad = sc.AnnData(X=X, obs=bc.set_index("barcode"),
                    var=feats.set_index("symbol"))
    ad.var_names_make_unique()
    pos = pos.set_index("barcode")
    ad.obs = ad.obs.join(pos[["in_tissue", "pxl_row_in_fullres",
                              "pxl_col_in_fullres"]])
    ad = ad[ad.obs["in_tissue"] == 1].copy()
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    coords = ad.obs[["pxl_col_in_fullres",
                     "pxl_row_in_fullres"]].to_numpy(float)
    return ad, coords


# ------------------------------------------------------- Moran's I (vectorised)
def morans_all(Z, coords, k=6, n_perm=499, seed=42):
    """Z: n_spots x n_genes dense, already column-centred.
    Returns (I, p_perm) arrays, two-sided permutation test."""
    from sklearn.neighbors import NearestNeighbors
    n, g = Z.shape
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()
    W = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    W = W.maximum(W.T)
    W = sparse.csr_matrix(W / W.sum(axis=1))
    S0 = W.sum()
    denom = (Z * Z).sum(axis=0)
    valid = denom > 0
    I = np.full(g, np.nan)
    P = np.full(g, np.nan)
    Iv = np.zeros(valid.sum())
    Zv = Z[:, valid]
    dv = denom[valid]
    Iv = (n / S0) * (Zv * (W @ Zv)).sum(axis=0) / dv
    rng = np.random.default_rng(seed)
    ge = np.zeros(Zv.shape[1], dtype=int)
    for _ in range(n_perm):
        p = rng.permutation(n)
        Zp = Zv[p]
        Ip = (n / S0) * (Zp * (W @ Zp)).sum(axis=0) / dv
        ge += (np.abs(Ip) >= np.abs(Iv))
    I[valid] = Iv
    P[valid] = (ge + 1) / (n_perm + 1)
    return I, P


def run_sample(tag, prefix, pos_file, species, treatment, timepoint,
               gene_sets, union_genes):
    out_csv = OUT / f"svg_{tag}.csv"
    if out_csv.exists():
        print(f"[skip] {tag}")
        return
    t0 = time.time()
    ad, coords = load_sample(prefix, pos_file)
    if species == "rat":
        mapping = map_orthologs(union_genes, ad.var_names)
    else:
        vs = set(ad.var_names)
        mapping = {g: (g if g in vs else None) for g in union_genes}
    sym2set = {}
    for sname, genes in gene_sets.items():
        for g in genes:
            sym2set.setdefault(g, []).append(sname)
    present = [(g, m) for g, m in mapping.items() if m and m in set(ad.var_names)]
    genes_h = [g for g, _ in present]
    syms = [m for _, m in present]
    X = ad[:, syms].X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    Z = X - X.mean(axis=0)
    I, P = morans_all(Z, coords)
    df = pd.DataFrame({
        "gene": genes_h, "sample": tag, "species": species,
        "treatment": treatment, "timepoint": timepoint,
        "symbol_in_sample": syms, "morans_i": I, "p_perm": P,
        "significant_p005": P < 0.05, "n_spots": ad.n_obs,
    })
    df["gene_sets"] = ["|".join(sym2set.get(g, [])) for g in genes_h]
    df.to_csv(out_csv, index=False)
    n_sig = int((P < 0.05).sum())
    print(f"[done] {tag}: spots={ad.n_obs} genes={len(genes_h)} "
          f"sig={n_sig} ({time.time()-t0:.0f}s)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", choices=["rat", "human"], required=True)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=None)
    args = ap.parse_args()
    gene_sets, union_genes = build_gene_sets()
    if args.species == "rat":
        samples = RAT_SAMPLES[args.start:args.end]
        for smp in samples:
            tr, tp = RAT_META[smp]
            run_sample(smp, str(RAT_DIR / smp), "tissue_positions_list.csv.gz",
                       "rat", tr, tp, gene_sets, union_genes)
    else:
        for tag, (pfx, tr) in HUM_SAMPLES.items():
            run_sample(tag, str(HUM_DIR / pfx), "tissue_positions.csv.gz",
                       "human", tr, "human", gene_sets, union_genes)


if __name__ == "__main__":
    main()
