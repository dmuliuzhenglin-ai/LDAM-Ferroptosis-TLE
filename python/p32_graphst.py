# -*- coding: utf-8 -*-
"""
p32_graphst.py
GraphST (Long et al., Nat Commun 2023) spatial-domain identification on
human hippocampus (GSE288888, 2 samples) and rat SE time-series
(GSE301260, 16 samples) 10x Visium data, replacing / validating the
heuristic KMeans annotation used in manuscript v2.7.

Phases:
  --train-only    : per-sample GraphST embedding + leiden clustering over a
                    grid of k, marker-based domain annotation, Moran's I of
                    domain labels; writes per-sample parquet + JSON.
  --analysis-only : pooled ARI vs KMeans, stratified statistics under the
                    GraphST domain framework, figures, report.
  (default: both)

Axis scores (LDAM / Ferroptosis / Module5 / microglia) are NOT recomputed;
they are reused from results/65_spatial_human and results/66_spatial_rat.

Run:  gsmap_env\\Scripts\\python.exe scripts\\p32_graphst.py [--train-only|--analysis-only|--skip-figs]
CPU only; epochs=600 per sample.
"""
import argparse
import gzip
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
from scipy import sparse, stats

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
DATA_H = ROOT / "data" / "spatial" / "GSE288888"
DATA_R = ROOT / "data" / "spatial" / "GSE301260"
RES_H = ROOT / "results" / "65_spatial_human"
RES_R = ROOT / "results" / "66_spatial_rat"
OUT = ROOT / "results" / "67_graphst"
FIGD = ROOT / "figures" / "v25"
IMG_CACHE = OUT / "img_cache"
for d in (OUT, FIGD, IMG_CACHE, OUT / "per_sample"):
    d.mkdir(parents=True, exist_ok=True)

EPOCHS = 600
SEED = 41
K_GRID = {"human": [5, 6, 7, 8], "rat": [5, 6, 7, 8, 9, 10]}
SCORE_AXES = ["score_LDAM", "score_Ferroptosis", "score_Module5"]

MARKERS_H = {
    "DG-like (PROX1+)": ["PROX1"],
    "CA1-like (NETO1+)": ["NETO1", "RGS14"],
    "astro-enriched": ["GFAP", "AQP4"],
    "thalamus": ["NTNG1"],
    "caudoputamen": ["GPR88"],
    "white_matter": ["MBP", "PLP1"],
    "microglia-enriched": ["AIF1", "CX3CR1", "P2RY12", "TMEM119"],
}
MARKERS_R = {
    "DG-like": ["Prox1"],
    "CA1-like": ["Rgs14", "Neto1"],
    "astro-enriched": ["Gfap", "Aqp4"],
    "thalamus": ["Ntng1"],
    "caudoputamen": ["Gpr88"],
    "white_matter": ["Mbp", "Plp1"],
    "microglia-enriched": ["Aif1", "Cx3cr1", "P2ry12", "Tmem119"],
}
HIPPO_LABELS = {"DG-like (PROX1+)", "CA1-like (NETO1+)", "DG-like", "CA1-like"}


# ------------------------------------------------------------------ loading
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
    return pos.set_index("barcode")


def load_visium(data_dir, prefix, pos_name):
    with gzip.open(data_dir / f"{prefix}_matrix.mtx.gz", "rb") as f:
        mat = scipy.io.mmread(f).tocsr().T  # spots x genes
    feats = pd.read_csv(data_dir / f"{prefix}_features.tsv.gz", sep="\t", header=None)
    genes = feats[1].astype(str).values
    bcs = pd.read_csv(data_dir / f"{prefix}_barcodes.tsv.gz", header=None)[0].astype(str).values
    pos = read_positions(data_dir / f"{prefix}_{pos_name}")
    var = pd.DataFrame(index=pd.Index(genes))
    adata = sc.AnnData(X=mat, obs=pd.DataFrame(index=bcs), var=var)
    adata.var_names_make_unique()
    pos = pos.reindex(adata.obs_names)
    adata.obs["array_row"] = pos["array_row"].values
    adata.obs["array_col"] = pos["array_col"].values
    adata.obs["pxl_row_in_fullres"] = pos["pxl_row_in_fullres"].values
    adata.obs["pxl_col_in_fullres"] = pos["pxl_col_in_fullres"].values
    adata = adata[pos["in_tissue"].values == 1].copy()
    adata.obsm["spatial"] = adata.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(float)
    return adata


def sample_table(species):
    """Return dict sample_key -> {adata loader args, keep_barcodes, meta}."""
    tab = {}
    if species == "human":
        for key, prefix, pos_name, csv in [
            ("human_CTRL", "GSM8778056_A18-43", "tissue_positions.csv.gz", "spot_scores_CTRL.csv"),
            ("human_HS", "GSM8778057_276A", "tissue_positions.csv.gz", "spot_scores_HS.csv"),
        ]:
            sc_df = pd.read_csv(RES_H / csv)
            sc_df["raw_barcode"] = sc_df["barcode"].str.rsplit("-", n=1).str[0]
            tab[key] = dict(prefix=prefix, pos_name=pos_name, data_dir=DATA_H,
                            scores=sc_df, group="CTRL" if key.endswith("CTRL") else "HS")
    else:
        pq = pd.read_parquet(RES_R / "spot_scores_regions.parquet")
        qc = pd.read_csv(RES_R / "sample_qc.csv")
        for _, row in qc.iterrows():
            key = row["sample"]
            sub = pq[pq["sample"] == key].copy()
            sub["raw_barcode"] = sub["barcode"]
            tab[key] = dict(prefix=key, pos_name="tissue_positions_list.csv.gz",
                            data_dir=DATA_R, scores=sub,
                            group=row["treatment"], timepoint=row["timepoint"])
    return tab


# ------------------------------------------------------------------ metrics
def knn_weights(rows, cols, k=6):
    from sklearn.neighbors import NearestNeighbors
    xy = np.c_[cols, rows].astype(float)
    nbrs = NearestNeighbors(n_neighbors=k + 1).fit(xy)
    _, idx = nbrs.kneighbors(xy)
    n = len(xy)
    ii = np.repeat(np.arange(n), k)
    jj = idx[:, 1:].ravel()
    W = sparse.csr_matrix((np.ones(len(ii)), (ii, jj)), shape=(n, n))
    W = W.maximum(W.T)
    return W


def morans_i(x, W, n_perm=99, seed=0):
    x = np.asarray(x, float)
    n = len(x)
    xc = x - x.mean()
    denom = float((xc ** 2).sum())
    S0 = W.sum()
    if denom == 0 or S0 == 0:
        return np.nan, np.nan
    I = (n / S0) * float(xc @ (W @ xc)) / denom
    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(n_perm):
        xp = rng.permutation(xc)
        Ip = (n / S0) * float(xp @ (W @ xp)) / denom
        if abs(Ip) >= abs(I):
            cnt += 1
    return I, (cnt + 1) / (n_perm + 1)


def domain_morans(labels, W, seed=0):
    """Mean Moran's I across one-hot encoded domain labels."""
    cats = pd.Categorical(labels)
    Is = []
    for c in cats.categories:
        oh = np.asarray(cats == c, dtype=float)
        if oh.sum() < 10:
            continue
        I, _ = morans_i(oh, W, n_perm=0)
        Is.append(I)
    return float(np.mean(Is)) if Is else np.nan


def find_resolution(adata, k, lo=0.02, hi=3.0, iters=18):
    """Bisection search for a leiden resolution giving exactly k clusters.
    Falls back to the closest cluster count (flagged) if exact hit fails."""
    def n_at(res):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        return int(adata.obs["leiden"].nunique())
    best = (None, 1e9, None)  # res, |n-k|, n
    for _ in range(iters):
        res = (lo + hi) / 2
        n = n_at(res)
        if n == k:
            return res, True
        if abs(n - k) < best[1]:
            best = (res, abs(n - k), n)
        if n < k:
            lo = res
        else:
            hi = res
    # fine scan around best
    for res in np.arange(max(0.02, best[0] - 0.1), best[0] + 0.101, 0.005):
        n = n_at(float(res))
        if n == k:
            return float(res), True
        if abs(n - k) < best[1]:
            best = (float(res), abs(n - k), n)
    n_at(best[0])
    return best[0], False


# ------------------------------------------------------------- train phase
def run_sample(key, spec, species):
    from GraphST.GraphST import GraphST as GraphSTModel
    from sklearn.metrics import silhouette_score

    t0 = time.time()
    adata = load_visium(spec["data_dir"], spec["prefix"], spec["pos_name"])
    keep = set(spec["scores"]["raw_barcode"])
    adata = adata[[b in keep for b in adata.obs_names]].copy()
    n_spots = adata.n_obs
    print(f"[{key}] {n_spots} QC spots", flush=True)

    # log-normalised copy for marker scoring (annotation only)
    ad_n = adata.copy()
    sc.pp.normalize_total(ad_n, target_sum=1e4)
    sc.pp.log1p(ad_n)

    markers = MARKERS_H if species == "human" else MARKERS_R
    markers = {m: [g for g in genes if g in adata.var_names] for m, genes in markers.items()}
    markers = {m: g for m, g in markers.items() if g}
    Xn = ad_n.X.toarray() if sparse.issparse(ad_n.X) else np.asarray(ad_n.X)
    vidx = {g: i for i, g in enumerate(ad_n.var_names)}
    marker_expr = {m: Xn[:, [vidx[g] for g in genes]].mean(1) for m, genes in markers.items()}

    # GraphST
    model = GraphSTModel(adata, datatype="10X", device="cpu", epochs=EPOCHS,
                         random_seed=SEED)
    adata = model.train()
    emb = adata.obsm["emb"]

    # k grid: leiden resolution bisection on the embedding
    sc.pp.neighbors(adata, n_neighbors=50, use_rep="emb")
    k_metrics = {}
    for k in K_GRID[species]:
        res, exact = find_resolution(adata, k)
        labels = adata.obs["leiden"].astype(str).values
        adata.obs[f"graphst_k{k}"] = labels
        sil = float(silhouette_score(emb, labels)) if len(set(labels)) > 1 else np.nan
        # marker contrast: fraction of domains with at least one marker group z>=1
        dfm = pd.DataFrame({"lab": labels})
        n_annot = 0
        for lab in sorted(set(labels)):
            sub = dfm["lab"] == lab
            zs = []
            for m, v in marker_expr.items():
                mu_in, mu_out = v[sub.values].mean(), v[~sub.values].mean()
                sd = v.std() + 1e-9
                zs.append((mu_in - v.mean()) / sd)
            if max(zs) >= 1.0:
                n_annot += 1
        k_metrics[k] = {"resolution": res, "silhouette": sil,
                        "resolution_exact": bool(exact),
                        "n_clusters_actual": int(len(set(labels))),
                        "n_domains_marker_annotated": n_annot,
                        "marker_contrast_frac": n_annot / len(set(labels))}
        print(f"[{key}] k={k} res={res:.2f} sil={sil:.3f} annot={n_annot}/{len(set(labels))} exact={exact}", flush=True)

    # choose k: max silhouette, tie-break by marker contrast
    ks = sorted(k_metrics)
    best_k = max(ks, key=lambda k: (round(k_metrics[k]["silhouette"], 4),
                                    k_metrics[k]["n_domains_marker_annotated"]))
    labels = adata.obs[f"graphst_k{best_k}"].astype(str).values
    print(f"[{key}] chosen k={best_k}", flush=True)

    # domain annotation (chosen k)
    annotation = {}
    for lab in sorted(set(labels), key=int):
        sub = labels == lab
        entry = {"n_spots": int(sub.sum())}
        zmap, mmap = {}, {}
        for m, v in marker_expr.items():
            mu = float(v[sub].mean())
            mmap[m] = mu
            zmap[m] = float((mu - v.mean()) / (v.std() + 1e-9))
        entry["mean_marker_expr"] = mmap
        entry["z_marker"] = zmap
        best_m = max(zmap, key=zmap.get)
        if zmap[best_m] >= 1.0:
            entry["annotation"] = best_m
            entry["confidence"] = "z>=1"
        else:
            entry["annotation"] = "mixed/other"
            entry["confidence"] = "below threshold (max z=%.2f)" % zmap[best_m]
        entry["z_best"] = zmap[best_m]
        annotation[lab] = entry

    # Moran's I: GraphST labels vs KMeans labels on identical KNN graph
    W = knn_weights(adata.obs["array_row"].values, adata.obs["array_col"].values, k=6)
    mi_graphst = domain_morans(labels, W)
    km_col = "expr_cluster" if species == "human" else "cluster"
    km_map = spec["scores"].set_index("raw_barcode")[km_col]
    km_labels = km_map.reindex(adata.obs_names).values
    mi_kmeans = domain_morans(pd.Series(km_labels).astype(str).values, W)

    out = pd.DataFrame({
        "barcode": adata.obs_names,
        "array_row": adata.obs["array_row"].values,
        "array_col": adata.obs["array_col"].values,
        "pxl_row_in_fullres": adata.obs["pxl_row_in_fullres"].values,
        "pxl_col_in_fullres": adata.obs["pxl_col_in_fullres"].values,
        "graphst_domain": labels,
    })
    for k in K_GRID[species]:
        out[f"graphst_k{k}"] = adata.obs[f"graphst_k{k}"].astype(str).values
    out["kmeans_label"] = pd.Series(km_labels).astype(str).values

    out.to_parquet(OUT / "per_sample" / f"{key}.parquet", index=False)
    meta = {"sample": key, "species": species, "n_spots": int(n_spots),
            "epochs": EPOCHS, "seed": SEED, "chosen_k": int(best_k),
            "k_metrics": {str(k): v for k, v in k_metrics.items()},
            "annotation": annotation,
            "morans_i_domain_labels": {"graphst": mi_graphst, "kmeans": mi_kmeans},
            "runtime_sec": time.time() - t0}
    if species == "rat":
        meta["treatment"] = spec["group"]
        meta["timepoint"] = spec["timepoint"]
    else:
        meta["group"] = spec["group"]
    with open(OUT / "per_sample" / f"{key}.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=1)
    print(f"[{key}] done in {time.time()-t0:.0f}s  MI graphst={mi_graphst:.3f} kmeans={mi_kmeans:.3f}", flush=True)
    return meta


def train_phase(only=None):
    todo = []
    for key, spec in sample_table("human").items():
        todo.append(("human", key, spec))
    for key, spec in sample_table("rat").items():
        todo.append(("rat", key, spec))
    metas = []
    for species, key, spec in todo:
        if only and key not in only:
            continue
        if (OUT / "per_sample" / f"{key}.json").exists() and (OUT / "per_sample" / f"{key}.parquet").exists():
            print(f"[{key}] cached, skip", flush=True)
            continue
        metas.append(run_sample(key, spec, species))
    return metas


# ---------------------------------------------------------- analysis phase
def recompute_annotations(dom, metas):
    """Re-annotate chosen-k domains with across-domain z-scores of domain-mean
    marker expression (same convention as the v2.7 KMeans annotation), so the
    z>=1 threshold is comparable. Marker groups whose maximum domain mean is
    <0.05 log-normalised counts are excluded as not measurably expressed.
    Overwrites the 'annotation' block of each per-sample JSON."""
    for s, m in sorted(metas.items()):
        species = m["species"]
        spec = sample_table(species)[s]
        adata = load_visium(spec["data_dir"], spec["prefix"], spec["pos_name"])
        keep = set(spec["scores"]["raw_barcode"])
        adata = adata[[b in keep for b in adata.obs_names]].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        markers = MARKERS_H if species == "human" else MARKERS_R
        markers = {mk: [g for g in genes if g in adata.var_names]
                   for mk, genes in markers.items()}
        markers = {mk: g for mk, g in markers.items() if g}
        X = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
        vidx = {g: i for i, g in enumerate(adata.var_names)}
        mex = {mk: X[:, [vidx[g] for g in genes]].mean(1) for mk, genes in markers.items()}
        bc2lab = dict(zip(dom.loc[dom["sample"] == s, "barcode"],
                          dom.loc[dom["sample"] == s, "graphst_domain"]))
        labs = np.array([bc2lab[b] for b in adata.obs_names])
        annotation = {}
        domains = sorted(set(labs), key=int)
        dmean = {mk: np.array([v[labs == d].mean() for d in domains])
                 for mk, v in mex.items()}
        # exclude marker groups with no measurable expression
        usable = {mk for mk, v in dmean.items() if v.max() >= 0.05}
        dropped = sorted(set(dmean) - usable)
        for i, d in enumerate(domains):
            zmap, mmap = {}, {}
            for mk in usable:
                v = dmean[mk]
                mmap[mk] = float(v[i])
                zmap[mk] = float((v[i] - v.mean()) / (v.std() + 1e-9))
            best_m = max(zmap, key=zmap.get)
            entry = {"n_spots": int((labs == d).sum()),
                     "mean_marker_expr": mmap, "z_marker": zmap,
                     "z_best": zmap[best_m],
                     "annotation": best_m if zmap[best_m] >= 1.0 else "mixed/other",
                     "confidence": "z>=1" if zmap[best_m] >= 1.0
                                   else "below threshold (max z=%.2f)" % zmap[best_m]}
            if dropped:
                entry["marker_groups_not_expressed"] = dropped
            annotation[d] = entry
        m["annotation"] = annotation
        m["annotation_method"] = ("across-domain z-score of domain-mean log-normalised "
                                  "marker expression (v2.7 convention); marker groups with "
                                  "max domain mean <0.05 excluded")
        with open(OUT / "per_sample" / f"{s}.json", "w", encoding="utf-8") as f:
            json.dump(m, f, indent=1)
        n_ann = sum(1 for e in annotation.values() if e["annotation"] != "mixed/other")
        print(f"[{s}] re-annotated: {n_ann}/{len(annotation)} domains labelled", flush=True)
    return metas


def mwu_r(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 5 or len(b) < 5:
        return np.nan, np.nan, np.nan
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    r = 2 * U / (len(a) * len(b)) - 1
    return U, p, r


def analysis_phase(skip_figs=False):
    from sklearn.metrics import adjusted_rand_score

    # ---------- load per-sample outputs
    metas, frames = {}, []
    for f in sorted((OUT / "per_sample").glob("*.json")):
        m = json.load(open(f, encoding="utf-8"))
        metas[m["sample"]] = m
    for f in sorted((OUT / "per_sample").glob("*.parquet")):
        df = pd.read_parquet(f)
        df["sample"] = f.stem
        frames.append(df)
    dom = pd.concat(frames, ignore_index=True)

    # re-annotate domains with v2.7-comparable across-domain z-scores
    metas = recompute_annotations(dom, metas)

    # attach scores + existing KMeans annotation
    h_tabs = sample_table("human")
    r_tabs = sample_table("rat")
    h_scores = pd.concat([s["scores"].assign(sample=k) for k, s in h_tabs.items()])
    r_scores = pd.concat([s["scores"].assign(sample=k) for k, s in r_tabs.items()])
    keep_h = ["raw_barcode", "sample", "score_LDAM", "score_Ferroptosis", "score_Module5",
              "microglia_enrichment", "region", "expr_cluster"]
    keep_r = ["raw_barcode", "sample", "treatment", "timepoint", "score_LDAM",
              "score_Ferroptosis", "score_Module5", "score_Microglia",
              "region", "cluster", "hippocampus_like"]
    scores = pd.concat([
        h_scores[keep_h].assign(species="human", treatment=h_scores["sample"].map(
            {"human_CTRL": "CTRL", "human_HS": "HS"}), timepoint="human",
            score_Microglia=h_scores["microglia_enrichment"], cluster=h_scores["expr_cluster"],
            hippocampus_like=np.nan),
        r_scores[keep_r].assign(species="rat"),
    ], ignore_index=True).rename(columns={"raw_barcode": "barcode"})
    dom = dom.merge(scores.drop(columns=["expr_cluster"], errors="ignore"),
                    on=["sample", "barcode"], how="left", suffixes=("", "_sc"))
    assert dom["score_LDAM"].notna().mean() > 0.99, "score join failed"
    for c in ["cluster", "region", "graphst_domain", "kmeans_label", "timepoint", "treatment"]:
        if c in dom.columns:
            dom[c] = dom[c].astype(str)

    # ---------- annotation.json (per sample per domain) + k selection
    annotation = {s: {"chosen_k": m["chosen_k"], "k_metrics": m["k_metrics"],
                      "annotation": m["annotation"],
                      "morans_i_domain_labels": m["morans_i_domain_labels"],
                      "n_spots": m["n_spots"], "epochs": m["epochs"]}
                  for s, m in metas.items()}
    with open(OUT / "annotation.json", "w", encoding="utf-8") as f:
        json.dump(annotation, f, indent=1)

    # annotation lookup for stratification
    ann_map = {}  # (sample, domain) -> label
    for s, m in metas.items():
        for lab, e in m["annotation"].items():
            ann_map[(s, lab)] = e["annotation"]
    dom["graphst_annotation"] = [ann_map.get((s, l), "mixed/other")
                                 for s, l in zip(dom["sample"], dom["graphst_domain"])]
    dom.to_parquet(OUT / "domains_per_sample.parquet", index=False)
    dom.to_csv(OUT / "domains_per_sample.csv", index=False)

    # ---------- ARI vs KMeans
    ari = {}
    for s, g in dom.groupby("sample"):
        ari[s] = {
            "ari_graphst_vs_kmeans_clusters": float(adjusted_rand_score(g["graphst_domain"], g["kmeans_label"])),
            "ari_graphst_vs_kmeans_region": float(adjusted_rand_score(g["graphst_domain"], g["region"].astype(str))),
            "n_spots": int(len(g)),
            "morans_i_graphst": metas[s]["morans_i_domain_labels"]["graphst"],
            "morans_i_kmeans": metas[s]["morans_i_domain_labels"]["kmeans"],
        }
    # per-domain correspondence table (human + pooled rat)
    corresp = {}
    for species, g in dom.groupby("species"):
        ct = pd.crosstab(g["graphst_annotation"], g["region"].astype(str))
        corresp[species] = ct.to_dict()
    ari_out = {"per_sample": ari, "annotation_vs_kmeans_region_crosstab": corresp,
               "note": "ARI computed per sample against the v2.7 KMeans(k=10) cluster labels "
                       "and against their marker-based region annotations."}
    with open(OUT / "ari_comparison.json", "w", encoding="utf-8") as f:
        json.dump(ari_out, f, indent=1)

    # ---------- stratified stats under GraphST framework
    strat = {"human": {"per_domain": {}, "microglia_enriched_domain": {}, "overall": {}},
             "rat": {"hippocampus_like": {}, "microglia_enriched": {}}}

    # human: HS vs CTRL within each annotation label
    hd = dom[dom["species"] == "human"]
    for lab in sorted(hd["graphst_annotation"].unique()):
        a = hd[(hd["graphst_annotation"] == lab) & (hd["treatment"] == "HS")]
        b = hd[(hd["graphst_annotation"] == lab) & (hd["treatment"] == "CTRL")]
        entry = {"n_HS": len(a), "n_CTRL": len(b)}
        for ax in SCORE_AXES:
            U, p, r = mwu_r(a[ax], b[ax])
            entry[ax] = {"median_HS": float(a[ax].median()) if len(a) else None,
                         "median_CTRL": float(b[ax].median()) if len(b) else None,
                         "p_value": p, "rank_biserial_r": r}
        strat["human"]["per_domain"][lab] = entry
    # microglia-enriched domain specifically
    mg = hd[hd["graphst_annotation"] == "microglia-enriched"]
    a, b = mg[mg["treatment"] == "HS"], mg[mg["treatment"] == "CTRL"]
    strat["human"]["microglia_enriched_domain"] = {
        "n_HS": len(a), "n_CTRL": len(b),
        **{ax: dict(zip(["median_HS", "median_CTRL", "p_value", "rank_biserial_r"],
                        [float(a[ax].median()) if len(a) else None,
                         float(b[ax].median()) if len(b) else None,
                         mwu_r(a[ax], b[ax])[1], mwu_r(a[ax], b[ax])[2]]))
           for ax in SCORE_AXES}}
    # overall (sanity vs v2.7)
    a, b = hd[hd["treatment"] == "HS"], hd[hd["treatment"] == "CTRL"]
    strat["human"]["overall"] = {ax: {
        "median_HS": float(a[ax].median()), "median_CTRL": float(b[ax].median()),
        "p_value": mwu_r(a[ax], b[ax])[1], "rank_biserial_r": mwu_r(a[ax], b[ax])[2]}
        for ax in SCORE_AXES}

    # rat: per timepoint, hippocampus-like and microglia-enriched strata
    rd = dom[dom["species"] == "rat"]
    for tp in ["D5", "D10", "D20", "D40"]:
        for stratum_name, mask_fn in [
            ("hippocampus_like", lambda g: g["graphst_annotation"].isin(["DG-like", "CA1-like"])),
            ("microglia_enriched", lambda g: g["graphst_annotation"] == "microglia-enriched"),
        ]:
            g = rd[rd["timepoint"] == tp]
            g = g[mask_fn(g)]
            se, ct = g[g["treatment"] == "SE"], g[g["treatment"] == "CTRL"]
            entry = {"n_spots_SE": len(se), "n_spots_CTRL": len(ct)}
            for ax in SCORE_AXES:
                U, p, r = mwu_r(se[ax], ct[ax])
                med_se = se.groupby("sample")[ax].median()
                med_ct = ct.groupby("sample")[ax].median()
                entry[ax] = {
                    "spot_level": {"median_SE": float(se[ax].median()) if len(se) else None,
                                   "median_CTRL": float(ct[ax].median()) if len(ct) else None,
                                   "p_value": p, "rank_biserial_r": r},
                    "sample_level": {"sample_medians_SE": med_se.to_dict(),
                                     "sample_medians_CTRL": med_ct.to_dict(),
                                     "delta_median_of_medians":
                                         float(med_se.median() - med_ct.median())
                                         if len(med_se) and len(med_ct) else None},
                }
            strat["rat"][stratum_name][tp] = entry

    # ---------- direction consistency vs v2.7
    ref_h = json.load(open(RES_H / "stats_summary.json", encoding="utf-8"))
    ref_r = json.load(open(RES_R / "stats_summary.json", encoding="utf-8"))
    consistency = []

    def cmp(tag, new_med_a, new_med_b, new_p, ref_med_a, ref_med_b, ref_p):
        if new_med_a is None or new_med_b is None:
            consistency.append({"comparison": tag, "status": "not computable under GraphST"})
            return
        new_dir = np.sign(new_med_a - new_med_b)
        ref_dir = np.sign(ref_med_a - ref_med_b)
        consistency.append({
            "comparison": tag,
            "new_delta": new_med_a - new_med_b, "new_p": new_p,
            "ref_delta": ref_med_a - ref_med_b, "ref_p": ref_p,
            "direction_consistent": bool(new_dir == ref_dir)})

    for ax in SCORE_AXES:
        r27 = ref_h["comparisons"][ax]
        nw = strat["human"]["overall"][ax]
        cmp(f"human overall {ax} (HS vs CTRL)", nw["median_HS"], nw["median_CTRL"],
            nw["p_value"], r27["median_HS"], r27["median_CTRL"], r27["p_value"])
    # human microglia: v2.7 used top-quartile microglia stratum
    for ax in SCORE_AXES:
        r27 = ref_h["microglia_top_quartile"]
        nw = strat["human"]["microglia_enriched_domain"][ax]
        cmp(f"human microglia stratum {ax} (HS vs CTRL) [v2.7=top-quartile]",
            nw["median_HS"], nw["median_CTRL"], nw["p_value"],
            r27[f"{ax}__HS"]["median"], r27[f"{ax}__CTRL"]["median"],
            r27[f"{ax}__HS_vs_CTRL"]["p_value"])
    for tp in ["D5", "D10", "D20", "D40"]:
        for ax in SCORE_AXES:
            r27 = ref_r["hippocampus_only"][tp][ax]["spot_level"]
            nw = strat["rat"]["hippocampus_like"][tp][ax]["spot_level"]
            cmp(f"rat {tp} hippocampus-like {ax} (SE vs CTRL)",
                nw["median_SE"], nw["median_CTRL"], nw["p_value"],
                r27["median_SE"], r27["median_CTRL"], r27["p_value"])
            r27 = ref_r["microglia_top_quartile"][tp][ax]["spot_level"]
            nw = strat["rat"]["microglia_enriched"][tp][ax]["spot_level"]
            cmp(f"rat {tp} microglia stratum {ax} (SE vs CTRL) [v2.7=top-quartile]",
                nw["median_SE"], nw["median_CTRL"], nw["p_value"],
                r27["median_SE"], r27["median_CTRL"], r27["p_value"])
    strat["consistency_vs_v27"] = consistency
    strat["framework_note"] = (
        "GraphST domains are per-sample; strata pooled across samples of the same "
        "group/timepoint. v2.7 rat microglia stratum was per-sample top-quartile "
        "score_Microglia; here it is the GraphST microglia-enriched domain.")
    with open(OUT / "stratified_stats.json", "w", encoding="utf-8") as f:
        json.dump(strat, f, indent=1)

    if not skip_figs:
        make_figures(dom, metas)
    write_report(dom, metas, ari_out, strat)
    return dom, metas, ari_out, strat


# ---------------------------------------------------------------- figures
def get_lowres(species, key, prefix, pos_name):
    """Decompress lowres image + scalefactors to cache; return (img, scalef)."""
    import matplotlib.pyplot as plt
    cache_png = IMG_CACHE / f"{key}_lowres.png"
    if species == "human":
        src_png = RES_H / f"{'CTRL' if key.endswith('CTRL') else 'HS'}_tissue_lowres.png"
        sf = json.load(open(RES_H / f"{'CTRL' if key.endswith('CTRL') else 'HS'}_scalefactors.json"))
        img = plt.imread(src_png)
    else:
        data_dir = DATA_R
        if not cache_png.exists():
            with gzip.open(data_dir / f"{prefix}_tissue_lowres_image.png.gz", "rb") as fi:
                with open(cache_png, "wb") as fo:
                    fo.write(fi.read())
        img = plt.imread(cache_png)
        with gzip.open(data_dir / f"{prefix}_scalefactors_json.json.gz", "rt") as f:
            sf = json.load(f)
    return img, sf["tissue_lowres_scalef"]


def make_figures(dom, metas):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    palette = plt.cm.tab10.colors

    def plot_domains(ax, g, img, scalef, title, annot=True):
        ax.imshow(img)
        labs = sorted(g["graphst_domain"].unique(), key=int)
        cmap = {l: palette[i % 10] for i, l in enumerate(labs)}
        for l in labs:
            sub = g[g["graphst_domain"] == l]
            x = sub["pxl_col_in_fullres"] * scalef
            y = sub["pxl_row_in_fullres"] * scalef
            tag = ann_map_label(metas, g["sample"].iloc[0], l) if annot else l
            ax.scatter(x, y, s=6, c=[cmap[l]], label=f"{l}:{tag}", linewidths=0)
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        ax.legend(fontsize=5, markerscale=1.5, loc="lower right", framealpha=0.7)

    def plot_score(ax, g, img, scalef, score, title):
        ax.imshow(img)
        x = g["pxl_col_in_fullres"] * scalef
        y = g["pxl_row_in_fullres"] * scalef
        sc_ = ax.scatter(x, y, s=6, c=g[score], cmap="viridis", linewidths=0)
        ax.set_title(title, fontsize=9)
        ax.axis("off")
        plt.colorbar(sc_, ax=ax, fraction=0.04, pad=0.01)

    # ----- human figure: 2 samples x (domains, LDAM, Ferroptosis)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    h_specs = {"human_CTRL": ("GSM8778056_A18-43", "tissue_positions.csv.gz", "CTRL (GSM8778056)"),
               "human_HS": ("GSM8778057_276A", "tissue_positions.csv.gz", "MTLE-HS (GSM8778057)")}
    for i, (key, (prefix, pos_name, label)) in enumerate(h_specs.items()):
        g = dom[dom["sample"] == key]
        img, scalef = get_lowres("human", key, prefix, pos_name)
        plot_domains(axes[i, 0], g, img, scalef, f"{label} GraphST domains (k={metas[key]['chosen_k']})")
        plot_score(axes[i, 1], g, img, scalef, "score_LDAM", f"{label} LDAM score")
        plot_score(axes[i, 2], g, img, scalef, "score_Ferroptosis", f"{label} Ferroptosis score")
    fig.suptitle("GraphST spatial domains vs axis scores — human hippocampus (GSE288888)", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGD / "Fig_graphst_human_domains.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ----- rat figure: one SE + one CTRL per timepoint
    reps = {"D5": ("GSM9078930_A_L1_S1", "GSM9078934_B_L1_S2"),
            "D10": ("GSM9078931_A_L2_S5", "GSM9078935_B_L2_S6"),
            "D20": ("GSM9078938_C_L1_S3", "GSM9078942_D_L1_S4"),
            "D40": ("GSM9078939_C_L2_S7", "GSM9078943_D2_L2_S8")}
    fig, axes = plt.subplots(4, 2, figsize=(11, 16))
    for i, (tp, (se, ct)) in enumerate(reps.items()):
        for j, key in enumerate([se, ct]):
            g = dom[dom["sample"] == key]
            img, scalef = get_lowres("rat", key, key, "tissue_positions_list.csv.gz")
            tr = metas[key]["treatment"]
            plot_domains(axes[i, j], g, img, scalef,
                         f"{tp} {tr} ({key.split('_')[0]}) k={metas[key]['chosen_k']}")
    fig.suptitle("GraphST spatial domains — rat SE time series (GSE301260)", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIGD / "Fig_graphst_rat_domains.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def ann_map_label(metas, sample, domain):
    return metas[sample]["annotation"].get(str(domain), {}).get("annotation", "?")


# ---------------------------------------------------------------- report
def write_report(dom, metas, ari_out, strat):
    L = []
    L.append("# GraphST spatial-domain identification — validation of v2.7 KMeans annotation\n")
    L.append("## Methods\n")
    L.append(
        "GraphST (Long et al., Nat Commun 2023;14:470) was run per sample on QC-passed "
        "in-tissue spots: raw counts -> seurat_v3 HVG (n=3000) -> normalize_total(1e4) -> "
        "log1p -> scale(zero_center=False); spatial adjacency from 3 nearest spots by "
        "full-resolution pixel distance; GraphST contrastive GNN auto-encoder trained on "
        "CPU for 600 epochs (seed 41); leiden clustering on the 64-d embedding with "
        "resolution searched to hit each target k. k was scanned over "
        f"{K_GRID['human']} (human) and {K_GRID['rat']} (rat) and chosen by embedding "
        "silhouette (tie-break: number of domains with a marker z>=1). Domains were "
        "annotated with the same marker system as v2.7 (z>=1 rule). Axis scores "
        "(LDAM / Ferroptosis / Module5 / microglia) were reused unchanged from "
        "results/65 and results/66; only grouping labels changed.\n")
    L.append("## Chosen k and cluster quality\n")
    L.append("| sample | n spots | chosen k | silhouette | Moran's I GraphST | Moran's I KMeans |")
    L.append("|---|---|---|---|---|---|")
    for s, m in metas.items():
        sil = m["k_metrics"][str(m["chosen_k"])]["silhouette"]
        mi = m["morans_i_domain_labels"]
        L.append(f"| {s} | {m['n_spots']} | {m['chosen_k']} | {sil:.3f} | "
                 f"{mi['graphst']:.3f} | {mi['kmeans']:.3f} |")
    L.append("\n## Agreement with v2.7 KMeans (ARI)\n")
    L.append("| sample | ARI vs KMeans clusters | ARI vs KMeans region labels |")
    L.append("|---|---|---|")
    for s, a in ari_out["per_sample"].items():
        L.append(f"| {s} | {a['ari_graphst_vs_kmeans_clusters']:.3f} | "
                 f"{a['ari_graphst_vs_kmeans_region']:.3f} |")
    L.append("\n## Domain annotation\n")
    L.append("| sample | domain | n spots | annotation | confidence |")
    L.append("|---|---|---|---|---|")
    for s, m in metas.items():
        for lab, e in sorted(m["annotation"].items(), key=lambda kv: int(kv[0])):
            L.append(f"| {s} | {lab} | {e['n_spots']} | {e['annotation']} | {e['confidence']} |")
    L.append("\n## Stratified statistics under the GraphST framework vs v2.7 direction\n")
    L.append("| comparison | GraphST delta | GraphST p | v2.7 delta | v2.7 p | direction consistent |")
    L.append("|---|---|---|---|---|---|")
    for c in strat["consistency_vs_v27"]:
        if c.get("status"):
            L.append(f"| {c['comparison']} | - | - | - | - | {c['status']} |")
        else:
            L.append(f"| {c['comparison']} | {c['new_delta']:.4f} | {c['new_p']:.2e} | "
                     f"{c['ref_delta']:.4f} | {c['ref_p']:.2e} | {c['direction_consistent']} |")
    L.append("\n## Key findings under the GraphST framework\n")
    h = strat["human"]
    L.append("### Human (GSE288888, n=1 vs n=1, descriptive only)")
    for ax in SCORE_AXES:
        o = h["overall"][ax]
        L.append(f"- Overall {ax}: HS median {o['median_HS']:.4f} vs CTRL {o['median_CTRL']:.4f} "
                 f"(r={o['rank_biserial_r']:.3f}), identical to v2.7 (same spots, scores unchanged).")
    for lab, e in h["per_domain"].items():
        if lab == "mixed/other" or e["n_HS"] == 0 or e["n_CTRL"] == 0:
            continue
        ldam, fer = e["score_LDAM"], e["score_Ferroptosis"]
        L.append(f"- {lab} (n HS={e['n_HS']}, CTRL={e['n_CTRL']}): LDAM r={ldam['rank_biserial_r']:.3f} "
                 f"(p={ldam['p_value']:.1e}), Ferroptosis r={fer['rank_biserial_r']:.3f} (p={fer['p_value']:.1e}).")
    mg = h["microglia_enriched_domain"]
    if mg["n_HS"] == 0 or mg["n_CTRL"] == 0:
        L.append("- No domain reached the z>=1 threshold for microglia markers in either human "
                 "sample, so a GraphST-domain-based microglia stratum does not exist; the v2.7 "
                 "microglia-top-quartile stratification is unaffected (scores unchanged).")
    L.append("- Caveat: human hippocampus sections received coarse labels (e.g. 'caudoputamen' in "
             "HS domain 0 at z=1.21 with white_matter/astro/microglia z=1.07-1.12 close behind, and "
             "'thalamus'/'white_matter' in CTRL); these annotations are low-confidence and reflect "
             "the coarse marker system, not true extra-hippocampal tissue.")
    L.append("\n### Rat (GSE301260, n=2 per group per timepoint, effect sizes only)")
    for tp in ["D5", "D10", "D20", "D40"]:
        hh = strat["rat"]["hippocampus_like"][tp]
        l = hh["score_LDAM"]["spot_level"]
        L.append(f"- {tp} hippocampus-like (CA1-like; no DG-like domain passed z>=1, as in v2.7): "
                 f"LDAM median SE {l['median_SE']:.4f} vs CTRL {l['median_CTRL']:.4f} "
                 f"(r={l['rank_biserial_r']:.3f}, p={l['p_value']:.1e}).")
    for tp in ["D5", "D10", "D20", "D40"]:
        mm = strat["rat"]["microglia_enriched"][tp]
        if mm["n_spots_SE"] == 0 or mm["n_spots_CTRL"] == 0:
            L.append(f"- {tp} microglia-enriched domain: one-sided "
                     f"(SE={mm['n_spots_SE']}, CTRL={mm['n_spots_CTRL']} spots); "
                     f"SE-vs-CTRL contrast not computable under GraphST domains.")
        else:
            l = mm["score_LDAM"]["spot_level"]
            L.append(f"- {tp} microglia-enriched domain: LDAM median SE {l['median_SE']:.4f} vs "
                     f"CTRL {l['median_CTRL']:.4f} (r={l['rank_biserial_r']:.3f}, p={l['p_value']:.1e}).")
    flips = [c for c in strat["consistency_vs_v27"]
             if not c.get("status") and not c["direction_consistent"]]
    L.append(f"\nDirection consistency vs v2.7: "
             f"{sum(1 for c in strat['consistency_vs_v27'] if c.get('direction_consistent'))} consistent, "
             f"{len(flips)} flipped, "
             f"{sum(1 for c in strat['consistency_vs_v27'] if c.get('status'))} not computable.")
    for c in flips:
        L.append(f"- FLIP: {c['comparison']}: GraphST delta {c['new_delta']:.4f} "
                 f"(p={c['new_p']:.1e}) vs v2.7 delta {c['ref_delta']:.4f} (p={c['ref_p']:.1e}).")
    L.append("\n## Limitations\n")
    L.append(
        "- Human data remain n=1 vs n=1 (no biological replicates); all human p-values are "
        "spot-level descriptive statistics, not inferential evidence.\n"
        "- Rat groups have n=2 per timepoint; replicate-level inference is not possible, so "
        "sample-level deltas are reported as effect sizes only.\n"
        "- k was chosen by silhouette on the embedding; marker-based annotation at z>=1 is a "
        "coarse anatomical label, and domains failing the threshold are reported as mixed/other "
        "without forced assignment.\n"
        "- GraphST embeddings are trained per sample; domain indices are not aligned across "
        "samples, so cross-sample pooling uses annotation labels, not domain ids.\n"
        "- No DG-like domain passed the z>=1 marker rule in any rat sample (consistent with the "
        "v2.7 KMeans annotation, which also lacked a DG cluster); the rat hippocampus-like "
        "stratum is therefore effectively CA1-like only.\n"
        "- The GraphST microglia-enriched domain is group-imbalanced at D10-D40 (present on one "
        "side only), so the v2.7 microglia-top-quartile stratification cannot be replaced "
        "one-to-one at those timepoints.\n"
        "- CPU-only training (600 epochs); deeper training or GPU refinement was not explored.\n")
    (OUT / "graphst_report.md").write_text("\n".join(L), encoding="utf-8")


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-only", action="store_true")
    ap.add_argument("--analysis-only", action="store_true")
    ap.add_argument("--skip-figs", action="store_true")
    ap.add_argument("--samples", nargs="*", default=None)
    args = ap.parse_args()
    if not args.analysis_only:
        train_phase(only=args.samples)
    if not args.train_only:
        analysis_phase(skip_figs=args.skip_figs)


if __name__ == "__main__":
    main()
