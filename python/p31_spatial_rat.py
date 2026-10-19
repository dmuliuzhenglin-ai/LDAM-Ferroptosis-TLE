# -*- coding: utf-8 -*-
"""
p31_spatial_rat.py
Rat epileptogenesis time-series spatial transcriptomics (10x Visium) pipeline.
Dataset: GSE301260 (Dufour A et al., Acta Neuropathol Commun. 2026;14:38).
Lithium-pilocarpine SE model; D5/D10/D20/D40 x SE/CTRL x 2 replicates.

Pipeline (parameterised: scans the data directory, analyses every sample whose
matrix.mtx.gz passes a full gzip integrity test; corrupt/missing samples are
recorded and skipped):

  1. Sample inventory + QC table (SOFT parsing, per-sample filtering).
  2. Human -> rat ortholog mapping of gene sets (title-case first, unmatched
     genes via mygene.info v3 API; cached to ortholog_map.json).
  3. Gene sets: LDAM signature, ferroptosis signature, Module5 (28 genes),
     microglia markers, hub single genes.
  4. Per-sample scanpy.tl.score_genes (ctrl_size=50).
  5. Coarse anatomical region annotation (pooled PCA + leiden/KMeans,
     marker-based cluster annotation; confidence reported).
  6. Statistics: spot-level Mann-Whitney + rank-biserial r per timepoint,
     sample-level median summaries (n<=2 per group -> effect sizes only,
     descriptive p where computable), hippocampus-only and
     microglia-top-quartile stratifications, KNN(k=6) Moran's I.
  7. Deliverables under results/66_spatial_rat/ and figures/v25/.

Run:  gsmap_env\\Scripts\\python.exe scripts\\p31_spatial_rat.py [--figs-only|--skip-figs]
"""
import argparse
import gzip
import json
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "spatial" / "GSE301260"
SOFT = ROOT / "data" / "spatial" / "GSE301260.soft.txt"
OUT = ROOT / "results" / "66_spatial_rat"
FIGD = ROOT / "figures" / "v25"
IMG_CACHE = OUT / "img_cache"
ORTHO_CACHE = OUT / "ortholog_map.json"
for d in (OUT, FIGD, IMG_CACHE):
    d.mkdir(parents=True, exist_ok=True)

MICROGLIA_MARKERS = ["AIF1", "CX3CR1", "P2RY12", "TMEM119", "C1QA", "CSF1R"]
HUB_GENES = ["LGALS3", "SPP1", "APOE", "FTH1"]
# Coarse anatomy after Dufour et al.: cortex / white matter / thalamus /
# hippocampus (plexiform, CA1, DG, CA2/3) / caudoputamen
REGION_MARKERS = {
    "DG-like": ["Prox1"],
    "CA1-like": ["Rgs14"],
    "thalamus": ["Ntng1"],
    "caudoputamen": ["Gpr88"],
    "white_matter": ["Mbp", "Plp1"],
}
AUX_MARKERS = ["Gfap", "Aqp4"]
SCORE_AXES = ["score_LDAM", "score_Ferroptosis", "score_Module5"]


# ---------------------------------------------------------------- SOFT parse
def parse_soft(path):
    """GSM -> {treatment, stage, slide, capture_area, title}.
    NOTE: in this SOFT file !Sample_title precedes !Sample_geo_accession,
    so titles are buffered and attached to the next accession."""
    samples, cur, pending_title = {}, None, None
    for line in Path(path).read_text(errors="replace").splitlines():
        line = line.strip()
        if line.startswith("!Sample_title"):
            pending_title = line.split("=")[1].strip()
            if cur and "title" not in samples.get(cur, {}):
                samples[cur]["title"] = pending_title
                pending_title = None
        elif line.startswith("!Sample_geo_accession"):
            cur = line.split("=")[1].strip()
            samples[cur] = {}
            if pending_title is not None:
                samples[cur]["title"] = pending_title
                pending_title = None
        elif cur and line.startswith("!Sample_characteristics"):
            val = line.split("=", 1)[1].strip()
            if ":" in val:
                k, v = val.split(":", 1)
                samples[cur][k.strip().lower().replace(" ", "_")] = v.strip()
    for gsm, d in samples.items():
        d["timepoint"] = f"D{d.get('stage', '?')}"
    return samples


def tp_sort(tps):
    """Natural time order D5 < D10 < D20 < D40."""
    return sorted(tps, key=lambda s: int(str(s).lstrip("D")))


# ------------------------------------------------------------------ loading
def gzip_ok(path):
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except Exception:
        return False


def read_positions(path):
    """Handle headerless (v1, tissue_positions_list) and headered (v2) files."""
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


def load_sample(prefix):
    with gzip.open(f"{prefix}_matrix.mtx.gz", "rb") as f:
        mat = scipy.io.mmread(f)
    X = sparse.csr_matrix(mat).T.tocsr()  # spots x genes
    feats = pd.read_csv(f"{prefix}_features.tsv.gz", sep="\t", header=None,
                        names=["id", "symbol", "ftype"])
    bc = pd.read_csv(f"{prefix}_barcodes.tsv.gz", header=None,
                     names=["barcode"])
    pos = read_positions(f"{prefix}_tissue_positions_list.csv.gz")
    ad = sc.AnnData(X=X, obs=bc.set_index("barcode"),
                    var=feats.set_index("symbol"))
    ad.var_names_make_unique()
    pos = pos.set_index("barcode")
    ad.obs = ad.obs.join(pos[["in_tissue", "array_row", "array_col",
                              "pxl_row_in_fullres", "pxl_col_in_fullres"]])
    n_raw = ad.n_obs
    ad = ad[ad.obs["in_tissue"] == 1].copy()
    ad.obsm["spatial"] = ad.obs[["pxl_col_in_fullres",
                                 "pxl_row_in_fullres"]].to_numpy(float)
    return ad, n_raw


def qc_normalize(ad):
    ad.var["mt"] = ad.var_names.str.match(r"(?i)^mt[-_]")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True)
    n_tissue = ad.n_obs
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    return ad, n_tissue


# ----------------------------------------------------------- ortholog mapping
def map_orthologs(human_genes, rat_symbols):
    """human symbol -> rat symbol. exact -> title-case -> case-insensitive ->
    mygene.info (Homo sapiens -> Rattus norvegicus via pantherdb ortholog).
    Cached on disk. Returns (mapping dict, provenance dict, network_ok)."""
    rat_lower = {}
    for s in rat_symbols:
        rat_lower.setdefault(s.lower(), s)
    cache = {}
    if ORTHO_CACHE.exists():
        cache = json.loads(ORTHO_CACHE.read_text(encoding="utf-8"))
    prov_sidecar = ORTHO_CACHE.with_name("ortholog_provenance.json")
    prov_cache = {}
    if prov_sidecar.exists():
        prov_cache = json.loads(prov_sidecar.read_text(encoding="utf-8"))
    mapping, provenance = {}, {}
    unmatched = []
    for g in human_genes:
        if g in cache:
            mapping[g] = cache[g]
            provenance[g] = prov_cache.get(g, "cache")
            continue
        if g in rat_symbols:
            mapping[g] = g
            provenance[g] = "exact"
        elif g.title() in rat_symbols:
            mapping[g] = g.title()
            provenance[g] = "title_case"
        elif g.lower() in rat_lower:
            mapping[g] = rat_lower[g.lower()]
            provenance[g] = "case_insensitive"
        else:
            unmatched.append(g)
    network_ok = None
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
                          "species": "human"},
                    timeout=30)
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
                              "scopes": "entrezgene",
                              "fields": "symbol", "species": "rat"},
                        timeout=30)
                    rr.raise_for_status()
                    ez2sym = {h.get("query"): h.get("symbol")
                              for h in rr.json() if "symbol" in h}
                    for g, ez in rat_entrez.items():
                        sym = ez2sym.get(ez)
                        if sym and sym in rat_symbols:
                            found[g] = sym
                        elif sym and sym.lower() in rat_lower:
                            found[g] = rat_lower[sym.lower()]
                time.sleep(0.34)  # mygene.info rate limit courtesy
            network_ok = True
            for g in unmatched:
                if g in found:
                    mapping[g] = found[g]
                    provenance[g] = "mygene_ortholog"
                else:
                    mapping[g] = None
                    provenance[g] = "no_rat_ortholog_found"
        except Exception as e:  # network failure -> title-case only
            network_ok = str(e)
            for g in unmatched:
                mapping[g] = None
                provenance[g] = "unmapped_network_failed"
    # persist cache (only confident mappings + explicit negatives)
    merged = dict(cache)
    for g, m in mapping.items():
        if g not in merged or merged[g] is None:
            merged[g] = m
    ORTHO_CACHE.write_text(json.dumps(merged, indent=1), encoding="utf-8")
    prov_merged = dict(prov_cache)
    for g, p in provenance.items():
        prov_merged.setdefault(g, p)
    prov_sidecar.write_text(json.dumps(prov_merged, indent=1), encoding="utf-8")
    return mapping, provenance, network_ok


# ------------------------------------------------------------------ statistics
def rank_biserial(x, y):
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    rbc = 2 * u / (len(x) * len(y)) - 1
    return float(u), float(p), float(rbc)


def knn_morans_i(x, coords, k=6, n_perm=499, seed=42):
    from sklearn.neighbors import NearestNeighbors
    n = len(x)
    if n < k + 5 or np.std(x) == 0:
        return None, None
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()
    W = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    W = W.maximum(W.T)
    W = sparse.csr_matrix(W / W.sum(axis=1))
    S0 = W.sum()
    z = x - x.mean()
    denom = z @ z
    I = (n / S0) * (z @ W @ z) / denom
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        zp = rng.permutation(z)
        Ip = (n / S0) * (zp @ W @ zp) / (zp @ zp)
        if abs(Ip) >= abs(I):
            ge += 1
    return float(I), float((ge + 1) / (n_perm + 1))


def compare(df, group_col, axes, label):
    """Spot-level MW (SE vs CTRL) + sample-level median summary."""
    out = {}
    se = df[df[group_col] == "SE"]
    ct = df[df[group_col] == "CTRL"]
    for ax in axes:
        x, y = se[ax].to_numpy(float), ct[ax].to_numpy(float)
        if len(x) == 0 or len(y) == 0:
            out[ax] = {"error": "empty group"}
            continue
        u, p, rbc = rank_biserial(x, y)
        med_se = se.groupby("sample", observed=True)[ax].median()
        med_ct = ct.groupby("sample", observed=True)[ax].median()
        out[ax] = {
            "spot_level": {
                "n_spots_SE": int(len(x)), "n_spots_CTRL": int(len(y)),
                "median_SE": float(np.median(x)), "median_CTRL": float(np.median(y)),
                "mannwhitney_U": u, "p_value": p, "rank_biserial_r": rbc},
            "sample_level": {
                "n_samples_SE": int(len(med_se)), "n_samples_CTRL": int(len(med_ct)),
                "sample_medians_SE": {k: float(v) for k, v in med_se.items()},
                "sample_medians_CTRL": {k: float(v) for k, v in med_ct.items()},
                "delta_median_of_medians": float(med_se.mean() - med_ct.mean()),
                "note": ("n<=2 per group; effect size only, no inferential "
                         "test at replicate level")}}
    return out


# ---------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figs-only", action="store_true")
    ap.add_argument("--skip-figs", action="store_true")
    args = ap.parse_args()
    if not args.figs_only:
        run_analysis()
    if not args.skip_figs:
        run_figures()


def run_analysis():
    soft = parse_soft(SOFT)
    print(f"[soft] {len(soft)} GSM entries parsed")

    # ---- inventory: scan directory, validate gzip ----
    mtx_files = sorted(DATA.glob("*_matrix.mtx.gz"))
    valid, excluded = {}, []
    for m in mtx_files:
        prefix = m.name.replace("_matrix.mtx.gz", "")
        gsm = prefix.split("_")[0]
        trio_ok = all((DATA / f"{prefix}_{s}").exists()
                      for s in ("barcodes.tsv.gz", "features.tsv.gz",
                                "tissue_positions_list.csv.gz"))
        if not trio_ok:
            excluded.append({"sample": prefix, "gsm": gsm,
                             "reason": "missing 10X trio file(s)"})
            continue
        if not gzip_ok(m):
            excluded.append({"sample": prefix, "gsm": gsm,
                             "reason": "matrix.mtx.gz failed gzip integrity test"})
            continue
        valid[gsm] = prefix
    print(f"[scan] valid={len(valid)} excluded={len(excluded)}")
    for e in excluded:
        print("  excluded:", e)

    # ---- load + QC ----
    ads, qc_rows = {}, []
    for gsm, prefix in sorted(valid.items()):
        meta = soft.get(gsm, {})
        ad, n_raw = load_sample(DATA / prefix)
        ad, n_tissue = qc_normalize(ad)
        ad.obs["sample"] = prefix
        ad.obs["gsm"] = gsm
        ad.obs["treatment"] = meta.get("treatment", "?")
        ad.obs["timepoint"] = meta.get("timepoint", "?")
        ad.obs_names = [f"{prefix}|{b}" for b in ad.obs_names]
        ads[prefix] = ad
        qc_rows.append({
            "gsm": gsm, "sample": prefix,
            "treatment": meta.get("treatment", "?"),
            "timepoint": meta.get("timepoint", "?"),
            "slide": meta.get("slide", "?"),
            "capture_area": meta.get("capture_area", "?"),
            "n_spots_raw": n_raw, "n_spots_in_tissue": n_tissue,
            "n_spots_after_qc": int(ad.n_obs),
            "n_genes_detected": int(ad.n_vars),
            "median_genes_per_spot": float(ad.obs["n_genes_by_counts"].median()),
            "median_counts_per_spot": float(ad.obs["total_counts"].median()),
            "median_mt_pct": float(ad.obs["pct_counts_mt"].median())})
        print(f"[load] {prefix} {meta.get('treatment')} {meta.get('timepoint')}: "
              f"{ad.n_obs} spots, {ad.n_vars} genes")
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(OUT / "sample_qc.csv", index=False)

    # ---- gene sets (human symbols) ----
    def read_gs(p):
        return [g.strip() for g in Path(p).read_text().splitlines() if g.strip()]
    ldam = read_gs(ROOT / "genesets" / "LDAM_geneset.txt")
    ferro = read_gs(ROOT / "genesets" / "Ferroptosis_geneset.txt")
    mod = pd.read_csv(ROOT / "results" / "59_metacell_wgcna" / "metacell_modules.csv")
    mod5 = mod.loc[mod["module"] == 5, "gene"].tolist()
    genesets_h = {"LDAM": ldam, "Ferroptosis": ferro, "Module5": mod5,
                  "Microglia": MICROGLIA_MARKERS, "HubGenes": HUB_GENES}
    geneset_sources = {
        "LDAM": "genesets/LDAM_geneset.txt (project LDAM signature)",
        "Ferroptosis": "genesets/Ferroptosis_geneset.txt (project ferroptosis set)",
        "Module5": "results/59_metacell_wgcna/metacell_modules.csv module==5",
        "Microglia": "canonical microglia markers (AIF1,CX3CR1,P2RY12,TMEM119,C1QA,CSF1R)",
        "HubGenes": "hub genes LGALS3,SPP1,APOE,FTH1"}

    # ---- ortholog mapping against union of rat symbols ----
    rat_symbols = set()
    for ad in ads.values():
        rat_symbols |= set(ad.var_names)
    all_h = sorted({g for gs in genesets_h.values() for g in gs} |
                   set(REGION_MARKERS) |
                   {m for v in REGION_MARKERS.values() for m in v} |
                   set(AUX_MARKERS))
    # REGION_MARKERS keys are already rat-case symbols; include as-is
    mapping, provenance, network = map_orthologs(all_h, rat_symbols)
    coverage, genesets_rat = {}, {}
    for name, gs in genesets_h.items():
        mapped = [mapping[g] for g in gs if mapping.get(g)]
        present = [g for g in mapped if g in rat_symbols]
        genesets_rat[name] = present
        coverage[name] = {
            "n_human_genes": len(gs),
            "n_mapped_to_rat": len(mapped),
            "n_present_in_matrix": len(present),
            "coverage_vs_matrix": round(len(present) / len(gs), 4),
            "unmapped": [g for g in gs if not mapping.get(g)]}
        print(f"[ortho] {name}: {len(present)}/{len(gs)} present")

    # ---- per-sample scoring ----
    for prefix, ad in ads.items():
        for name in ["LDAM", "Ferroptosis", "Module5", "Microglia"]:
            gs = [g for g in genesets_rat[name] if g in ad.var_names]
            if len(gs) >= 5:
                sc.tl.score_genes(ad, gs, ctrl_size=50,
                                  score_name=f"score_{name}", random_state=0)
            else:
                ad.obs[f"score_{name}"] = np.nan
                print(f"[warn] {prefix} {name}: only {len(gs)} genes, skip score")
        for g in HUB_GENES:
            rg = mapping.get(g)
            if rg and rg in ad.var_names:
                v = ad[:, rg].X
                ad.obs[f"expr_{g}"] = (v.toarray().ravel()
                                       if sparse.issparse(v) else np.asarray(v).ravel())
            else:
                ad.obs[f"expr_{g}"] = np.nan
    print("[score] done")

    # ---- region annotation (pooled clustering + marker z-score) ----
    combined = sc.concat(list(ads.values()), label="batch", join="inner",
                         index_unique=None)
    print(f"[merge] {combined.n_obs} spots, {combined.n_vars} shared genes")
    sc.pp.highly_variable_genes(combined, n_top_genes=2000, flavor="seurat",
                                batch_key="sample")
    clu = combined[:, combined.var["highly_variable"]].copy()
    sc.pp.scale(clu, max_value=10)
    sc.pp.pca(clu, n_comps=30, random_state=0)
    sc.pp.neighbors(clu, n_neighbors=15, n_pcs=30, random_state=0)
    cluster_method = "leiden"
    try:
        sc.tl.leiden(clu, resolution=0.8, random_state=0)
    except Exception as e:
        from sklearn.cluster import KMeans
        km = KMeans(n_clusters=10, n_init=10, random_state=0)
        clu.obs["leiden"] = km.fit_predict(clu.obsm["X_pca"]).astype(str)
        cluster_method = f"KMeans(k=10) fallback: {type(e).__name__}"
    combined.obs["cluster"] = clu.obs["leiden"].values
    print(f"[cluster] method={cluster_method}, "
          f"n_clusters={combined.obs['cluster'].nunique()}")

    # marker expression per cluster (log-normalised), z-score across clusters
    marker_genes = {lab: [g for g in gs if g in combined.var_names]
                    for lab, gs in REGION_MARKERS.items()}
    cl_means = {}
    for cl in sorted(combined.obs["cluster"].unique(),
                     key=lambda x: int(x)):
        m = (combined.obs["cluster"] == cl).to_numpy()
        vals = {}
        for lab, gs in marker_genes.items():
            if gs:
                sub = combined[m, :][:, gs].X
                if sparse.issparse(sub):
                    sub = sub.toarray()
                vals[lab] = float(np.mean(sub))
            else:
                vals[lab] = 0.0
        cl_means[cl] = vals
    zmat = pd.DataFrame(cl_means).T
    zmat = (zmat - zmat.mean()) / zmat.std().replace(0, np.nan)
    cl_annot, region_rows = {}, []
    for cl in zmat.index:
        zr = zmat.loc[cl]
        best, zmax = zr.idxmax(), zr.max()
        if pd.notna(zmax) and zmax >= 1.0:
            annot, conf = best, float(zmax)
        else:
            annot, conf = "mixed/other", float(zmax) if pd.notna(zmax) else 0.0
        cl_annot[cl] = annot
        n = int((combined.obs["cluster"] == cl).sum())
        region_rows.append({"cluster": cl, "annotation": annot,
                            "n_spots": n, "z_best": conf,
                            **{f"mean_{k}": float(zmat.loc[cl, k])
                               for k in zmat.columns},
                            **{f"z_{k}": float(zr[k]) for k in zmat.columns}})
    combined.obs["region"] = combined.obs["cluster"].map(cl_annot)
    combined.obs["region_z_best"] = combined.obs["cluster"].map(
        {r["cluster"]: r["z_best"] for r in region_rows})
    combined.obs["hippocampus_like"] = combined.obs["region"].isin(
        ["DG-like", "CA1-like"])
    region_df = pd.DataFrame(region_rows)
    region_df.to_csv(OUT / "region_annotation.csv", index=False)
    print("[region] " + "; ".join(
        f"c{r['cluster']}={r['annotation']}(n={r['n_spots']},z={r['z_best']:.2f})"
        for r in region_rows))

    # ---- spot table ----
    obs = combined.obs.copy()
    keep_cols = (["gsm", "sample", "treatment", "timepoint", "array_row",
                  "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres",
                  "n_genes_by_counts", "total_counts", "pct_counts_mt",
                  "cluster", "region", "region_z_best", "hippocampus_like"] +
                 SCORE_AXES + ["score_Microglia"] +
                 [f"expr_{g}" for g in HUB_GENES])
    spot_df = obs[keep_cols].copy()
    spot_df.insert(0, "barcode", [i.split("|", 1)[1] for i in spot_df.index])
    spot_df.to_parquet(OUT / "spot_scores_regions.parquet", index=False)
    print(f"[save] spot_scores_regions.parquet: {len(spot_df)} spots")

    # ---- statistics ----
    stats_out = {
        "dataset": "GSE301260 (Dufour et al. 2026) rat lithium-pilocarpine SE Visium",
        "pipeline": "scripts/p31_spatial_rat.py",
        "n_samples_valid": len(valid),
        "excluded_samples": excluded,
        "geneset_sources": geneset_sources,
        "ortholog_coverage": coverage,
        "ortholog_provenance_counts": {
            p: sum(1 for v in provenance.values() if v == p)
            for p in sorted(set(provenance.values()))},
        "mygene_network_status": ("ok" if network is True else
                                  (f"failed: {network}" if network else "not_needed")),
        "cluster_method": cluster_method,
        "scoring": ("per-sample scanpy.tl.score_genes ctrl_size=50 on "
                    "log-normalised counts (target_sum=1e4)"),
        "per_timepoint": {}, "hippocampus_only": {},
        "microglia_top_quartile": {}, "morans_i": {}, "timecourse": {},
        "region_annotation": region_rows,
    }
    axes4 = SCORE_AXES + ["score_Microglia"]
    tps = tp_sort(spot_df["timepoint"].unique())
    for tp in tps:
        sub = spot_df[spot_df["timepoint"] == tp]
        stats_out["per_timepoint"][tp] = {
            "n_samples": {t: int(sub.loc[sub["treatment"] == t, "sample"].nunique())
                          for t in ("SE", "CTRL")},
            "n_spots": {t: int((sub["treatment"] == t).sum())
                        for t in ("SE", "CTRL")},
            **compare(sub, "treatment", axes4, tp)}
        hc = sub[sub["hippocampus_like"]]
        stats_out["hippocampus_only"][tp] = {
            "n_spots": {t: int((hc["treatment"] == t).sum())
                        for t in ("SE", "CTRL")},
            **compare(hc, "treatment", SCORE_AXES, tp)}
        # microglia top quartile within each sample, then pooled
        mg_parts = []
        for s, sdf in sub.groupby("sample", observed=True):
            thr = sdf["score_Microglia"].quantile(0.75)
            mg_parts.append(sdf[sdf["score_Microglia"] >= thr])
        mg = pd.concat(mg_parts) if mg_parts else sub.iloc[0:0]
        stats_out["microglia_top_quartile"][tp] = {
            "definition": "score_Microglia >= per-sample 75th percentile",
            "n_spots": {t: int((mg["treatment"] == t).sum())
                        for t in ("SE", "CTRL")},
            **compare(mg, "treatment", SCORE_AXES, tp)}
        print(f"[stats] {tp} done")

    # ---- time-course trajectories (all spots / hippocampus-like / mg top-Q4) ----
    mg_mask = pd.Series(False, index=spot_df.index)
    for s, sdf in spot_df.groupby("sample", observed=True):
        thr = sdf["score_Microglia"].quantile(0.75)
        mg_mask.loc[sdf.index[sdf["score_Microglia"] >= thr]] = True
    layers = {"all_spots": spot_df,
              "hippocampus_like": spot_df[spot_df["hippocampus_like"]],
              "microglia_top_quartile": spot_df[mg_mask]}
    for layer_name, ldf in layers.items():
        tc = {}
        for ax in SCORE_AXES:
            ax_rec = {}
            for tr in ("SE", "CTRL"):
                pts = {}
                for tp in tps:
                    cell = ldf[(ldf["timepoint"] == tp) &
                               (ldf["treatment"] == tr)]
                    if len(cell) == 0:
                        continue
                    med = cell.groupby("sample", observed=True)[ax].median()
                    pts[tp] = {
                        "n_samples": int(len(med)),
                        "sample_medians": {k: float(v) for k, v in med.items()},
                        "mean_of_sample_medians": float(med.mean()),
                        "spot_median": float(cell[ax].median()),
                        "n_spots": int(len(cell))}
                ax_rec[tr] = pts
                # descriptive trend: delta between last and first timepoint
                ks = [k for k in tps if k in pts]
                if len(ks) >= 2:
                    d = (pts[ks[-1]]["mean_of_sample_medians"]
                         - pts[ks[0]]["mean_of_sample_medians"])
                    ax_rec[f"{tr}_delta_{ks[-1]}_vs_{ks[0]}"] = float(d)
            tc[ax] = ax_rec
        stats_out["timecourse"][layer_name] = tc
    print("[timecourse] done")

    moran_feats = axes4 + [f"expr_{g}" for g in HUB_GENES]
    for prefix, ad in ads.items():
        sub = spot_df[spot_df["sample"] == prefix]
        coords = sub[["pxl_col_in_fullres",
                      "pxl_row_in_fullres"]].to_numpy(float)
        for feat in moran_feats:
            x = sub[feat].to_numpy(float)
            if np.isnan(x).any():
                continue
            I, p = knn_morans_i(x, coords)
            if I is not None:
                stats_out["morans_i"][f"{feat}__{prefix}"] = {"I": I, "p_perm": p}
        print(f"[moran] {prefix} done")

    with open(OUT / "stats_summary.json", "w") as f:
        json.dump(stats_out, f, indent=2)
    print("[save] stats_summary.json")
    write_report(stats_out, qc, region_df)
    print("ANALYSIS DONE")


# -------------------------------------------------------------------- report
def write_report(st, qc, region_df):
    L = []
    A = L.append
    A("# Rat spatial transcriptomics (GSE301260) - LDAM x ferroptosis in epileptogenesis\n")
    A("## 1. Dataset and samples\n")
    A(f"- Dataset: {st['dataset']}")
    A(f"- Valid samples analysed: **{st['n_samples_valid']}** "
      f"(gzip-validated 10X matrices).")
    if st["excluded_samples"]:
        A(f"- **Excluded samples ({len(st['excluded_samples'])}):**")
        for e in st["excluded_samples"]:
            A(f"  - `{e['sample']}`: {e['reason']}")
    A("\nPer-timepoint sample sizes (biological replicates):")
    for tp, d in st["per_timepoint"].items():
        A(f"- {tp}: SE n={d['n_samples']['SE']}, CTRL n={d['n_samples']['CTRL']} "
          f"({d['n_spots']['SE']} vs {d['n_spots']['CTRL']} in-tissue QC-passed spots)")
    A("\n## 2. Gene sets and human->rat ortholog mapping\n")
    A("Mapping: exact rat symbol -> title-case -> case-insensitive -> "
      "mygene.info (pantherdb orthologs). Cached in `ortholog_map.json`. "
      f"Network status: **{st['mygene_network_status']}**.")
    A("\n| Gene set | Source | Human genes | Present in rat matrix | Coverage |")
    A("|---|---|---|---|---|")
    for name, cov in st["ortholog_coverage"].items():
        A(f"| {name} | {st['geneset_sources'][name]} | {cov['n_human_genes']} "
          f"| {cov['n_present_in_matrix']} | {cov['coverage_vs_matrix']:.1%} |")
    A("\nMapping provenance counts: "
      + ", ".join(f"{k}={v}" for k, v in st["ortholog_provenance_counts"].items()))
    A("\n## 3. Scoring and region annotation\n")
    A(f"- {st['scoring']}.")
    A(f"- Region annotation: pooled clustering ({st['cluster_method']}) on top-2000 "
      "HVGs; clusters annotated by z-scored marker means "
      "(Prox1=DG, Rgs14=CA1, Ntng1=thalamus, Gpr88=caudoputamen, Mbp/Plp1=white "
      "matter; z>=1 required, else 'mixed/other').")
    A("\n| Cluster | Annotation | n spots | best z |")
    A("|---|---|---|---|")
    for r in st["region_annotation"]:
        A(f"| {r['cluster']} | {r['annotation']} | {r['n_spots']} | {r['z_best']:.2f} |")
    A("\n## 4. Core statistics: SE vs CTRL per timepoint\n")
    for tp, d in st["per_timepoint"].items():
        A(f"\n### {tp} (spot-level Mann-Whitney; rank-biserial r>0 = higher in SE)\n")
        A("| Score | median SE | median CTRL | U | p | rank-biserial r |")
        A("|---|---|---|---|---|---|")
        for ax, res in d.items():
            if not isinstance(res, dict) or "spot_level" not in res:
                continue
            s = res["spot_level"]
            A(f"| {ax} | {s['median_SE']:.4f} | {s['median_CTRL']:.4f} "
              f"| {s['mannwhitney_U']:.3g} | {s['p_value']:.3g} | {s['rank_biserial_r']:+.4f} |")
        A("\nSample-level medians (no inferential test; n<=2 per group):")
        for ax, res in d.items():
            if not isinstance(res, dict) or "sample_level" not in res:
                continue
            s = res["sample_level"]
            A(f"- {ax}: SE medians { {k: round(v,4) for k,v in s['sample_medians_SE'].items()} }, "
              f"CTRL medians { {k: round(v,4) for k,v in s['sample_medians_CTRL'].items()} }, "
              f"delta(median-of-medians) = {s['delta_median_of_medians']:+.4f}")
    A("\n## 5. Stratified analyses\n")
    A("\n### Hippocampus-like spots only (DG-like + CA1-like clusters)")
    for tp, d in st["hippocampus_only"].items():
        A(f"\n**{tp}** (spots: SE={d['n_spots']['SE']}, CTRL={d['n_spots']['CTRL']}):")
        for ax in SCORE_AXES:
            res = d.get(ax, {})
            if "spot_level" in res:
                s = res["spot_level"]
                A(f"- {ax}: rank-biserial r = {s['rank_biserial_r']:+.4f}, "
                  f"p = {s['p_value']:.3g}")
            else:
                A(f"- {ax}: not estimable ({res.get('error','no spots')})")
    A("\n### Microglia-enriched spots (per-sample top quartile of microglia score)")
    for tp, d in st["microglia_top_quartile"].items():
        A(f"\n**{tp}** (spots: SE={d['n_spots']['SE']}, CTRL={d['n_spots']['CTRL']}):")
        for ax in SCORE_AXES:
            res = d.get(ax, {})
            if "spot_level" in res:
                s = res["spot_level"]
                A(f"- {ax}: rank-biserial r = {s['rank_biserial_r']:+.4f}, "
                  f"p = {s['p_value']:.3g}")
            else:
                A(f"- {ax}: not estimable ({res.get('error','no spots')})")
    A("\n## 6. Time-course trajectories (D5 -> D10 -> D20 -> D40)\n")
    A("Mean of per-sample median scores per timepoint (n=2 samples per "
      "group per timepoint; effect sizes descriptive).")
    layer_titles = {"all_spots": "All in-tissue spots",
                    "hippocampus_like": "Hippocampus-like spots (DG-like + CA1-like)",
                    "microglia_top_quartile": "Microglia-enriched top-quartile spots"}
    for layer, tc in st.get("timecourse", {}).items():
        A(f"\n### {layer_titles.get(layer, layer)}\n")
        A("| Score | Group | " + " | ".join(tp_sort(
            {tp for ax in tc.values() for tr in ('SE', 'CTRL')
             for tp in ax.get(tr, {})})) + " | delta last-first |")
        A("|---|---|" + "---|" * (len(tp_sort(
            {tp for ax in tc.values() for tr in ('SE', 'CTRL')
             for tp in ax.get(tr, {})})) + 1))
        tps_all = tp_sort({tp for ax in tc.values() for tr in ("SE", "CTRL")
                           for tp in ax.get(tr, {})})
        for ax, ax_rec in tc.items():
            for tr in ("SE", "CTRL"):
                pts = ax_rec.get(tr, {})
                cells = []
                for tp in tps_all:
                    if tp in pts:
                        cells.append(f"{pts[tp]['mean_of_sample_medians']:+.4f}")
                    else:
                        cells.append("-")
                dk = [k for k in ax_rec if k.startswith(f"{tr}_delta_")]
                dv = f"{ax_rec[dk[0]]:+.4f}" if dk else "-"
                A(f"| {ax.replace('score_','')} | {tr} | " + " | ".join(cells)
                  + f" | {dv} |")
    A("\n## 7. Spatial autocorrelation (Moran's I, KNN k=6, 499 permutations)\n")
    A("| Sample | Feature | I | perm p |")
    A("|---|---|---|---|")
    for k, v in st["morans_i"].items():
        feat, sample = k.split("__")
        A(f"| {sample} | {feat} | {v['I']:.4f} | {v['p_perm']:.3g} |")
    A("\n## 8. Limitations (read before interpreting)\n")
    A("- Two biological replicates per group per timepoint (2 SE + 2 CTRL, "
      "balanced). All spot-level p-values are descriptive: spots within a "
      "sample are spatially autocorrelated and not independent replicates. "
      "Replicate-level inference is impossible at n=2; only effect sizes "
      "(rank-biserial r, delta of sample medians) are reported.")
    A("- Region annotation is marker-heuristic (no squidpy/manual pathology "
      "review); 'mixed/other' clusters carry low confidence and hippocampus "
      "subfield calls (DG/CA1) should be treated as approximate.")
    A("- Human->rat ortholog coverage is <100% for large sets (see table); "
      "unmapped genes (e.g. human-specific lncRNAs such as AL163541.1) are "
      "listed in stats_summary.json.")
    A("- score_genes is run per sample; scores are internally normalised but "
      "cross-sample comparability relies on the shared control-gene-binning "
      "scheme, not on batch correction.")
    (OUT / "analysis_report.md").write_text("\n".join(L), encoding="utf-8")
    print("[save] analysis_report.md")


# ------------------------------------------------------------------- figures
def _load_img(gsm_prefix):
    """Decompress hires image + scalefactor into img_cache; return (img, sf)."""
    src = DATA / f"{gsm_prefix}_tissue_hires_image.png.gz"
    sf_src = DATA / f"{gsm_prefix}_scalefactors_json.json.gz"
    dst = IMG_CACHE / f"{gsm_prefix}_hires.png"
    try:
        if not dst.exists():
            with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        import matplotlib.image as mpimg
        img = mpimg.imread(dst)
        sf = json.load(gzip.open(sf_src, "rt"))["tissue_hires_scalef"]
        return img, sf
    except Exception as e:
        print(f"[fig] image unavailable for {gsm_prefix}: {e}")
        return None, None


def run_figures():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    spot = pd.read_parquet(OUT / "spot_scores_regions.parquet")
    st = json.loads((OUT / "stats_summary.json").read_text())

    # ---- Fig 1: spatial maps, representative SE vs CTRL pairs ----
    pairs = {}
    for tp in tp_sort(spot["timepoint"].unique()):
        sub = spot[spot["timepoint"] == tp]
        se_s = sorted(sub.loc[sub["treatment"] == "SE", "sample"].unique())
        ct_s = sorted(sub.loc[sub["treatment"] == "CTRL", "sample"].unique())
        if se_s and ct_s:
            pairs[tp] = (se_s[0], ct_s[0])
    tp_picks = [tp for tp in ("D10", "D40", "D20", "D5") if tp in pairs][:2]
    feats = [("score_LDAM", "LDAM score"), ("score_Ferroptosis", "Ferroptosis score"),
             ("score_Microglia", "Microglia score"), ("expr_LGALS3", "Lgals3 expr")]
    ncols = 2 * len(tp_picks)
    fig, axes = plt.subplots(4, ncols, figsize=(5.5 * ncols, 16),
                             squeeze=False)
    for row, (feat, lab) in enumerate(feats):
        for pi, tp in enumerate(tp_picks):
            se_id, ct_id = pairs[tp]
            ids = [se_id, ct_id]
            vals = spot.loc[spot["sample"].isin(ids), feat]
            vmin, vmax = np.nanpercentile(vals, [2, 98])
            for col, sid in enumerate(ids):
                ax = axes[row, 2 * pi + col]
                sub = spot[spot["sample"] == sid]
                img, sf = _load_img(sid)
                if img is not None:
                    ax.imshow(img)
                    x = sub["pxl_col_in_fullres"] * sf
                    y = sub["pxl_row_in_fullres"] * sf
                else:
                    x, y = sub["array_col"], -sub["array_row"]
                sca = ax.scatter(x, y, c=sub[feat], cmap="viridis", s=14,
                                 vmin=vmin, vmax=vmax, linewidths=0)
                tr = sub["treatment"].iloc[0]
                ax.set_title(f"{lab} | {sid.split('_')[0]} ({tr}, {tp})",
                             fontsize=9)
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_aspect("equal")
                fig.colorbar(sca, ax=ax, fraction=0.035, pad=0.01)
    fig.suptitle("GSE301260 rat Visium - representative SE vs CTRL pairs "
                 f"({', '.join(tp_picks)})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    f1 = FIGD / "Fig_spatial_rat_maps.png"
    fig.savefig(f1, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {f1}")

    # ---- Fig 2: per-timepoint SE vs CTRL distributions ----
    tps = tp_sort(spot["timepoint"].unique())
    fig, axes = plt.subplots(len(tps), 3, figsize=(13, 4.2 * len(tps)),
                             squeeze=False)
    for i, tp in enumerate(tps):
        sub = spot[spot["timepoint"] == tp]
        for j, ax_name in enumerate(SCORE_AXES):
            ax = axes[i, j]
            data = [sub.loc[sub["treatment"] == t, ax_name].dropna().to_numpy()
                    for t in ("SE", "CTRL")]
            parts = ax.violinplot(data, showmedians=True, widths=0.7)
            for pc, c in zip(parts["bodies"], ["#d1495b", "#2e86ab"]):
                pc.set_facecolor(c); pc.set_alpha(0.75)
            for k, d in enumerate(data):
                ax.scatter(np.full(len(d), k + 1) +
                           np.random.default_rng(0).normal(0, 0.04, len(d)),
                           d, s=1.5, c="k", alpha=0.25, linewidths=0)
            res = st["per_timepoint"][tp].get(ax_name, {})
            s = res.get("spot_level", {})
            ax.set_xticks([1, 2]); ax.set_xticklabels(["SE", "CTRL"])
            n_se = st["per_timepoint"][tp]["n_samples"]["SE"]
            n_ct = st["per_timepoint"][tp]["n_samples"]["CTRL"]
            ax.set_title(f"{tp} | {ax_name.replace('score_','')}\n"
                         f"r_rb={s.get('rank_biserial_r', float('nan')):+.3f}, "
                         f"p={s.get('p_value', float('nan')):.2e} "
                         f"(n={n_se}v{n_ct})", fontsize=9)
    fig.suptitle("SE vs CTRL spot-level score distributions "
                 "(rank-biserial r from Mann-Whitney; n = biological replicates)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    f2 = FIGD / "Fig_spatial_rat_stats.png"
    fig.savefig(f2, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {f2}")

    # ---- Fig 3: time-course trajectories (box + sample-median lines) ----
    layers = [("all_spots", "All spots", spot),
              ("hippocampus_like", "Hippocampus-like spots",
               spot[spot["hippocampus_like"]])]
    mg_mask = pd.Series(False, index=spot.index)
    for s, sdf in spot.groupby("sample", observed=True):
        thr = sdf["score_Microglia"].quantile(0.75)
        mg_mask.loc[sdf.index[sdf["score_Microglia"] >= thr]] = True
    layers.append(("microglia_top_quartile",
                   "Microglia top-quartile spots", spot[mg_mask]))
    colors = {"SE": "#d1495b", "CTRL": "#2e86ab"}
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(3, 3, figsize=(16, 13), squeeze=False)
    for i, (lkey, ltitle, ldf) in enumerate(layers):
        for j, ax_name in enumerate(SCORE_AXES):
            ax = axes[i, j]
            xs = np.arange(len(tps))
            for t_off, tr in [(-0.17, "SE"), (0.17, "CTRL")]:
                box_data, line_x, line_y = [], [], []
                for k, tp in enumerate(tps):
                    cell = ldf[(ldf["timepoint"] == tp) &
                               (ldf["treatment"] == tr)][ax_name].dropna()
                    box_data.append(cell.to_numpy())
                    med = (ldf[(ldf["timepoint"] == tp) &
                               (ldf["treatment"] == tr)]
                           .groupby("sample", observed=True)[ax_name]
                           .median())
                    if len(med):
                        line_x.append(k + t_off)
                        line_y.append(med.mean())
                        ax.scatter(np.full(len(med), k + t_off)
                                   + rng.normal(0, 0.02, len(med)),
                                   med.to_numpy(), s=28, zorder=5,
                                   facecolor="white", edgecolor=colors[tr],
                                   linewidth=1.4)
                bp = ax.boxplot(box_data, positions=xs + t_off, widths=0.28,
                                showfliers=False, patch_artist=True,
                                medianprops={"color": "k", "lw": 1},
                                whiskerprops={"color": colors[tr]},
                                capprops={"color": colors[tr]},
                                boxprops={"facecolor": colors[tr],
                                          "alpha": 0.45,
                                          "edgecolor": colors[tr]})
                ax.plot(line_x, line_y, color=colors[tr], lw=1.8,
                        marker="o", ms=4, label=tr if (i == 0 and j == 0)
                        else None, zorder=4)
            ax.set_xticks(xs); ax.set_xticklabels(tps)
            ax.set_title(f"{ax_name.replace('score_','')} | {ltitle}",
                         fontsize=10)
            ax.tick_params(labelsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", fontsize=10)
    fig.suptitle("Score trajectories across epileptogenesis (D5-D40): "
                 "box = spot-level distribution; line/dots = mean and "
                 "individual sample medians (n=2 per group)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    f3 = FIGD / "Fig_spatial_rat_timecourse.png"
    fig.savefig(f3, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {f3}")
    print("FIGURES DONE")


if __name__ == "__main__":
    main()
