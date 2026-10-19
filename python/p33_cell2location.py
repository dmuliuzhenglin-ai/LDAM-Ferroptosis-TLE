# -*- coding: utf-8 -*-
"""
p33_cell2location.py
cell2location spatial deconvolution of human hippocampus Visium (GSE288888)
using GSE190452 human temporal lobe snRNA-seq as reference.
Microglia are split into leiden substates; the LDAM-enriched substate
(leiden 35, highest mean LDAM_score) is mapped as cell type "MG_LDAMhi".
CPU only. Chunked training with checkpoints so the run can resume after
external kills. All numbers are written to results/68_cell2location/.
"""
import os, sys, json, time, gzip, warnings, argparse
warnings.filterwarnings("ignore")
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as st
import matplotlib
matplotlib.use("Agg")   # cell2location's filter_genes calls plt.show(); must be headless
import anndata as ad
import scanpy as sc
import scvi
from cell2location.models import RegressionModel, Cell2location
from cell2location.utils.filtering import filter_genes

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "results", "68_cell2location")
FIGDIR = os.path.join(BASE, "figures", "v25")
SPATIAL = os.path.join(BASE, "data", "spatial", "GSE288888")
S65 = os.path.join(BASE, "results", "65_spatial_human")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

scvi.settings.seed = 0
np.random.seed(0)

REF_H5AD = os.path.join(BASE, "results", "43_scDRS", "gse190452_annotated.h5ad")
MAX_CELLS_PER_LABEL = 1000
REF_MAX_EPOCHS = 250
REF_CHUNK = 50
GENE_CAP = 6000
N_CELLS_PER_LOCATION = 8
DETECTION_ALPHA = 200
MAP_TIME_BUDGET_S = 1500      # per-slice wall-clock budget for mapping training
MAP_MIN_EPOCHS = 1000
MAP_MAX_EPOCHS = 20000
MAP_CHUNK = 250
MAP_PROBE_EPOCHS = 10

SAMPLES = {
    "CTRL": "GSM8778056_A18-43",
    "HS":   "GSM8778057_276A",
}
LDAM_LABEL = "MG_LDAMhi"      # leiden 35 substate


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def save_json(obj, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, default=str)


def load_json(path, default=None):
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return default


# ------------------------------------------------------- chunked training
def train_chunked(mod, model_cls, adata, ckpt_dir, total_epochs, chunk,
                  train_kwargs, state, tag, init_kwargs=None):
    """Train in chunks, checkpointing after each chunk so a kill loses <chunk."""
    os.makedirs(ckpt_dir, exist_ok=True)
    state_f = os.path.join(ckpt_dir, "chunk_state.json")
    if state is None:
        state = {"epochs_done": 0, "history": []}
    done = state["epochs_done"]
    hist = state["history"]
    if done > 0:
        log(f"{tag}: resuming from checkpoint at epoch {done}")
        mod = model_cls.load(ckpt_dir, adata=adata)
    elif mod is None:
        mod = model_cls(adata, **(init_kwargs or {}))
    t0 = time.time()
    while done < total_epochs:
        step = min(chunk, total_epochs - done)
        mod.train(max_epochs=step, accelerator="cpu", **train_kwargs)
        h = mod.history["elbo_train"]["elbo_train"].values.tolist()
        hist.extend(h)
        done += step
        mod.save(ckpt_dir, overwrite=True)
        state.update({"epochs_done": done, "history": hist})
        save_json(state, state_f)
        log(f"{tag}: epoch {done}/{total_epochs} checkpointed, "
            f"ELBO={h[-1]:.5g}")
    return mod, time.time() - t0, np.array(hist)


# ---------------------------------------------------------------- reference
def prepare_reference():
    log("loading reference " + REF_H5AD)
    adata = ad.read_h5ad(REF_H5AD)
    adata.X = adata.layers["counts"].copy()
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)

    # --- two-level label: celltype, microglia split into leiden substates
    lab = adata.obs["celltype"].astype(str).copy()
    mg_mask = lab == "Microglia"
    mg_leiden = adata.obs.loc[mg_mask, "leiden"].astype(str)
    mg_lab = "MG_c" + mg_leiden
    mg_lab = mg_lab.where(mg_lab != "MG_c35", LDAM_LABEL)
    lab.loc[mg_mask] = mg_lab
    adata.obs["c2l_label"] = pd.Categorical(lab)
    log("labels: " + json.dumps(adata.obs["c2l_label"].value_counts().to_dict()))

    sub = adata.obs.loc[mg_mask].groupby("c2l_label", observed=True)["LDAM_score"].agg(
        ["count", "mean", "median"])
    sub.round(6).reset_index().to_json(
        os.path.join(OUT, "ref_microglia_substates.json"),
        orient="records", force_ascii=False)

    # --- stratified subsample
    keep = []
    rng = np.random.default_rng(0)
    for lb, idx in adata.obs.groupby("c2l_label", observed=True).indices.items():
        if len(idx) > MAX_CELLS_PER_LABEL:
            idx = np.sort(rng.choice(idx, MAX_CELLS_PER_LABEL, replace=False))
        keep.append(idx)
    keep = np.sort(np.concatenate(keep))
    adata = adata[keep].copy()
    log(f"subsampled reference: {adata.n_obs} cells")

    sel = filter_genes(adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03,
                       nonz_mean_cutoff=1.12)
    if hasattr(sel, "sum") and not isinstance(sel, pd.Index):
        log(f"filter_genes retained {int(sel.sum())} / {adata.n_vars}")
    else:
        sel = pd.Index(sel)
        log(f"filter_genes retained {len(sel)} / {adata.n_vars}")
    adata = adata[:, sel].copy()
    mt = adata.var_names.str.upper().str.startswith("MT-")
    adata = adata[:, ~mt].copy()

    spat_genes = set()
    for prefix in SAMPLES.values():
        with gzip.open(os.path.join(SPATIAL, prefix + "_features.tsv.gz"), "rt") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    spat_genes.add(parts[1])
    common = adata.var_names.intersection(pd.Index(sorted(spat_genes)))
    log(f"genes shared with Visium: {len(common)}")
    adata = adata[:, common].copy()
    if adata.n_vars > GENE_CAP:
        tot = np.asarray(adata.X.sum(axis=0)).ravel()
        top = np.argsort(-tot)[:GENE_CAP]
        adata = adata[:, np.sort(top)].copy()
    log(f"reference genes for training: {adata.n_vars}")

    RegressionModel.setup_anndata(adata=adata, batch_key="sample",
                                  labels_key="c2l_label")
    adata.write(os.path.join(OUT, "ref_adata_prepared.h5ad"))
    return adata


def run_reference():
    sig_csv = os.path.join(OUT, "reference_signatures.csv")
    if os.path.exists(sig_csv):
        log("loading existing reference signatures")
        return pd.read_csv(sig_csv, index_col=0)

    prep = os.path.join(OUT, "ref_adata_prepared.h5ad")
    ckpt = os.path.join(OUT, "ckpt_ref_model")
    state = load_json(os.path.join(ckpt, "chunk_state.json"))
    if state is not None and os.path.exists(prep):
        log("resuming: loading prepared reference anndata")
        adata = ad.read_h5ad(prep)
        RegressionModel.setup_anndata(adata=adata, batch_key="sample",
                                      labels_key="c2l_label")
    else:
        adata = prepare_reference()

    mod, train_s, hist = train_chunked(
        None, RegressionModel, adata, ckpt, REF_MAX_EPOCHS, REF_CHUNK,
        dict(batch_size=2500, lr=0.002), state, "reference")
    log(f"reference trained; ELBO first={hist[0]:.4g} last={hist[-1]:.4g} "
        f"(n_epochs={len(hist)})")

    adata = mod.export_posterior(
        adata, sample_kwargs={"num_samples": 1000, "batch_size": 2500})
    fn = adata.uns["mod"]["factor_names"]
    inf_aver = adata.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{i}" for i in fn]].copy()
    inf_aver.columns = fn
    inf_aver.to_csv(sig_csv)
    adata.write(os.path.join(OUT, "ref_adata_with_signatures.h5ad"))

    sub = pd.read_json(os.path.join(OUT, "ref_microglia_substates.json"))
    meta = {
        "ref_source": REF_H5AD,
        "ref_study": "GSE190452 (Chen et al. 2023, Nat Neurosci, PMID 36941428) human temporal lobe TLE snRNA-seq",
        "ref_cells_used": int(adata.n_obs),
        "ref_genes_used": int(adata.n_vars),
        "labels": adata.obs["c2l_label"].value_counts().to_dict(),
        "microglia_substate_LDAM_score": sub.to_dict("records"),
        "max_cells_per_label": MAX_CELLS_PER_LABEL,
        "ref_max_epochs": REF_MAX_EPOCHS,
        "ref_train_seconds_this_run": round(train_s, 1),
        "ref_elbo_first": float(hist[0]),
        "ref_elbo_last": float(hist[-1]),
        "ref_elbo_last10_mean": float(np.mean(hist[-10:])),
    }
    save_json(meta, os.path.join(OUT, "reference_meta.json"))
    log("reference stage done")
    return inf_aver


# ------------------------------------------------------------------ spatial
def clean_abund_columns(df):
    import re
    ren = {}
    for c in df.columns:
        m = re.match(r"^q(\d\d)[_q\d]*cell_abundance_w_sf_(.+)$", c)
        if m:
            ren[c] = f"q{m.group(1)}_{m.group(2)}"
    return df.rename(columns=ren)


def load_visium(prefix, tag):
    d = SPATIAL
    with gzip.open(os.path.join(d, prefix + "_matrix.mtx.gz"), "rb") as f:
        from scipy.io import mmread
        X = mmread(f).tocsr().T   # matrix is genes x barcodes -> spots x genes
    feats = pd.read_csv(gzip.open(os.path.join(d, prefix + "_features.tsv.gz"), "rt"),
                        sep="\t", header=None, names=["ensembl", "symbol", "ftype"])
    bcs = pd.read_csv(gzip.open(os.path.join(d, prefix + "_barcodes.tsv.gz"), "rt"),
                      header=None, names=["barcode"])["barcode"].values
    pos = pd.read_csv(gzip.open(os.path.join(d, prefix + "_tissue_positions.csv.gz"), "rt"))
    pos = pos.set_index("barcode").loc[bcs]
    in_tissue = pos["in_tissue"].values == 1
    X = X[in_tissue]
    pos = pos[in_tissue]
    bcs2 = bcs[in_tissue]

    var = feats[["symbol"]].copy()
    var.index = feats["symbol"].values
    var.index.name = None
    a = ad.AnnData(X=sp.csr_matrix(X.astype(np.float32)), var=var)
    a.var_names_make_unique()
    a.obs_names = [f"{bc}-{tag}" for bc in bcs2]
    a.obs["array_row"] = pos["array_row"].values
    a.obs["array_col"] = pos["array_col"].values
    a.obs["pxl_row_in_fullres"] = pos["pxl_row_in_fullres"].values
    a.obs["pxl_col_in_fullres"] = pos["pxl_col_in_fullres"].values
    a.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values
    a = a[:, ~a.var_names.str.upper().str.startswith("MT-")].copy()
    log(f"{tag}: {a.n_obs} in-tissue spots x {a.n_vars} genes")
    return a


def map_sample(tag, inf_aver):
    abund_f = os.path.join(OUT, f"abundance_{tag}.csv.gz")
    meta_f = os.path.join(OUT, f"mapping_meta_{tag}.json")
    if os.path.exists(abund_f) and os.path.exists(meta_f):
        log(f"{tag}: loading existing abundance")
        return clean_abund_columns(pd.read_csv(abund_f)), load_json(meta_f)

    vis = load_visium(SAMPLES[tag], tag)
    genes = inf_aver.index.intersection(vis.var_names)
    log(f"{tag}: {len(genes)} signature genes present in slice")
    vis = vis[:, genes].copy()
    sig = inf_aver.loc[genes].copy()
    Cell2location.setup_anndata(adata=vis)
    kwargs = dict(N_cells_per_location=N_CELLS_PER_LOCATION,
                  detection_alpha=DETECTION_ALPHA)

    ckpt = os.path.join(OUT, f"ckpt_map_{tag}")
    state = load_json(os.path.join(ckpt, "chunk_state.json"))
    if state is not None and state["epochs_done"] > 0:
        planned = state["planned_epochs"]
        mod = None
    else:
        state = {"epochs_done": 0, "history": []}
        # probe timing on a throwaway model
        probe = Cell2location(vis, cell_state_df=sig, **kwargs)
        t0 = time.time()
        probe.train(max_epochs=MAP_PROBE_EPOCHS, lr=0.002, accelerator="cpu")
        dt = (time.time() - t0) / MAP_PROBE_EPOCHS
        del probe
        planned = int(max(MAP_MIN_EPOCHS,
                          min(MAP_MAX_EPOCHS, MAP_TIME_BUDGET_S / max(dt, 1e-6))))
        state["planned_epochs"] = planned
        log(f"{tag}: probe {dt:.2f}s/epoch -> planned {planned} epochs")

    mod, train_s, hist = train_chunked(
        None, Cell2location, vis, ckpt, planned, MAP_CHUNK,
        dict(lr=0.002), state, f"map-{tag}", init_kwargs=dict(cell_state_df=sig, **kwargs))
    log(f"{tag}: mapping trained; ELBO first={hist[0]:.4g} last={hist[-1]:.4g}")

    vis = mod.export_posterior(
        vis, sample_kwargs={"batch_size": 2500},
        use_quantiles=True, add_to_obsm=["q05", "q50", "q95"])

    rows = []
    for q in ["q05", "q50", "q95"]:
        df = vis.obsm[f"{q}_cell_abundance_w_sf"].copy()
        df.columns = [f"{q}_{c.replace(f'{q}_cell_abundance_w_sf_', '')}"
                      for c in df.columns]
        rows.append(df)
    abund = pd.concat(rows, axis=1)
    abund.insert(0, "barcode", vis.obs_names)
    abund.insert(1, "sample", tag)
    abund.insert(2, "array_row", vis.obs["array_row"].values)
    abund.insert(3, "array_col", vis.obs["array_col"].values)
    abund.insert(4, "pxl_row_in_fullres", vis.obs["pxl_row_in_fullres"].values)
    abund.insert(5, "pxl_col_in_fullres", vis.obs["pxl_col_in_fullres"].values)
    abund.to_csv(abund_f, index=False, compression="gzip")
    vis.write(os.path.join(OUT, f"vis_{tag}_mapped.h5ad"))
    meta = {
        "n_spots": int(vis.n_obs),
        "n_signature_genes": int(len(genes)),
        "map_epochs": int(planned),
        "map_train_seconds_this_run": round(train_s, 1),
        "map_elbo_first": float(hist[0]),
        "map_elbo_last": float(hist[-1]),
        "map_elbo_last10_mean": float(np.mean(hist[-10:])),
        "N_cells_per_location": N_CELLS_PER_LOCATION,
        "detection_alpha": DETECTION_ALPHA,
    }
    save_json(meta, meta_f)
    return abund, meta


# --------------------------------------------------------------- downstream
def morans_i(x, rows, cols, n_perm=499, seed=0):
    pts = np.column_stack([rows, cols]).astype(float)
    d2 = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1)
    W = ((d2 > 0) & (d2 <= 2.05 ** 2)).astype(float)  # adjacent hex neighbours
    S0 = W.sum()
    n = len(x)
    z = x - x.mean()

    def mi(v):
        return (n / S0) * (W * np.outer(v, v)).sum() / (v * v).sum()

    I = mi(z)
    rng = np.random.default_rng(seed)
    perm = np.array([mi(z[rng.permutation(n)]) for _ in range(n_perm)])
    p = (1 + np.sum(np.abs(perm) >= abs(I))) / (n_perm + 1)
    return float(I), float(p), int(W.sum() / 2)


def rank_biserial(x_hs, x_ctrl):
    U, p = st.mannwhitneyu(x_hs, x_ctrl, alternative="two-sided")
    n1, n2 = len(x_hs), len(x_ctrl)
    r_hs_higher = 2 * U / (n1 * n2) - 1
    return float(U), float(p), float(r_hs_higher)


def downstream(abund_all):
    log("downstream analysis")
    scores = {}
    for tag in SAMPLES:
        s = pd.read_csv(os.path.join(S65, f"spot_scores_{tag}.csv"))
        scores[tag] = s.set_index("barcode")

    stats = {"samples": {}, "HS_vs_CTRL": {},
             "correlation_with_geneset_scores": {},
             "morans_I": {}, "hotspots_top_quartile": {}}
    merged = {}
    for tag in SAMPLES:
        a = abund_all[tag].set_index("barcode")
        overlap = [c for c in ["array_row", "array_col", "pxl_row_in_fullres",
                               "pxl_col_in_fullres"] if c in scores[tag].columns]
        j = a.join(scores[tag].drop(columns=overlap), how="inner")
        log(f"{tag}: {len(j)} spots joined with gene-set scores "
            f"(of {len(a)} mapped / {len(scores[tag])} scored)")
        merged[tag] = j
        col = f"q50_{LDAM_LABEL}"
        stats["samples"][tag] = {
            "n_spots_mapped": int(len(a)),
            "n_spots_joined": int(len(j)),
            f"q50_{LDAM_LABEL}_median": float(j[col].median()),
            f"q50_{LDAM_LABEL}_mean": float(j[col].mean()),
            f"q50_{LDAM_LABEL}_frac_nonzero": float((j[col] > 0).mean()),
        }
        for sc_col in ["score_LDAM", "score_Ferroptosis"]:
            rho, p = st.spearmanr(j[col], j[sc_col])
            stats["correlation_with_geneset_scores"].setdefault(sc_col, {})[tag] = {
                "spearman_rho": float(rho), "p_value": float(p), "n": int(len(j))}
        mg_cols = [c for c in a.columns if c.startswith("q50_MG")]
        j["q50_total_microglia"] = j[mg_cols].sum(axis=1)
        rho, p = st.spearmanr(j["q50_total_microglia"], j["microglia_enrichment"])
        stats["correlation_with_geneset_scores"].setdefault(
            "total_microglia_vs_microglia_enrichment", {})[tag] = {
            "spearman_rho": float(rho), "p_value": float(p), "n": int(len(j))}
        # relative (compositional) abundance of the LDAM-high substate
        ct_cols = [c for c in a.columns if c.startswith("q50_")
                   and not c.endswith(("array_row", "array_col"))]
        j[f"rel_{LDAM_LABEL}"] = j[col] / j[ct_cols].sum(axis=1)
        stats["samples"][tag][f"rel_{LDAM_LABEL}_median"] = float(
            j[f"rel_{LDAM_LABEL}"].median())
        stats["samples"][tag][f"rel_{LDAM_LABEL}_mean"] = float(
            j[f"rel_{LDAM_LABEL}"].mean())
        for sc_col in ["score_LDAM", "score_Ferroptosis"]:
            rho, p = st.spearmanr(j[f"rel_{LDAM_LABEL}"], j[sc_col])
            stats["correlation_with_geneset_scores"].setdefault(
                f"rel_{LDAM_LABEL}_vs_{sc_col}", {})[tag] = {
                "spearman_rho": float(rho), "p_value": float(p),
                "n": int(len(j))}
        I, pI, nedges = morans_i(j[col].values, j["array_row"].values,
                                 j["array_col"].values)
        stats["morans_I"][tag] = {"I": I, "perm_p": pI, "n_perm": 499,
                                  "n_neighbour_pairs": nedges}
        thr = j[col].quantile(0.75)
        hot = j[j[col] >= thr]
        reg = hot["region"].value_counts(normalize=True).round(4).to_dict() \
            if "region" in hot.columns else {}
        stats["hotspots_top_quartile"][tag] = {
            "threshold_q75": float(thr), "n_hotspots": int(len(hot)),
            "region_fraction_of_hotspots": reg,
            "hotspot_median_score_LDAM": float(hot["score_LDAM"].median()),
            "nonhotspot_median_score_LDAM": float(
                j.loc[j[col] < thr, "score_LDAM"].median()),
        }

    xh = merged["HS"][f"q50_{LDAM_LABEL}"].values
    xc = merged["CTRL"][f"q50_{LDAM_LABEL}"].values
    U, p, r = rank_biserial(xh, xc)
    stats["HS_vs_CTRL"][f"q50_{LDAM_LABEL}"] = {
        "median_HS": float(np.median(xh)), "median_CTRL": float(np.median(xc)),
        "mean_HS": float(np.mean(xh)), "mean_CTRL": float(np.mean(xc)),
        "mannwhitney_U": U, "p_value": p,
        "rank_biserial_r_positive_means_HS_higher": r,
        "note": "n=1 slice per group; descriptive proof-of-concept only, "
                "no biological replication",
    }
    mh = merged["HS"]["q50_total_microglia"].values
    mc = merged["CTRL"]["q50_total_microglia"].values
    U2, p2, r2 = rank_biserial(mh, mc)
    stats["HS_vs_CTRL"]["q50_total_microglia"] = {
        "median_HS": float(np.median(mh)), "median_CTRL": float(np.median(mc)),
        "mannwhitney_U": U2, "p_value": p2,
        "rank_biserial_r_positive_means_HS_higher": r2,
        "note": "n=1 slice per group; descriptive proof-of-concept only",
    }
    rel = f"rel_{LDAM_LABEL}"
    rh = merged["HS"][rel].values
    rc = merged["CTRL"][rel].values
    U3, p3, r3 = rank_biserial(rh, rc)
    stats["HS_vs_CTRL"][rel] = {
        "median_HS": float(np.median(rh)), "median_CTRL": float(np.median(rc)),
        "mean_HS": float(np.mean(rh)), "mean_CTRL": float(np.mean(rc)),
        "mannwhitney_U": U3, "p_value": p3,
        "rank_biserial_r_positive_means_HS_higher": r3,
        "note": "relative abundance = q50_MG_LDAMhi / sum of q50 over all "
                "mapped cell types per spot; n=1 vs n=1 descriptive",
    }

    save_json(stats, os.path.join(OUT, "stats_summary.json"))
    make_figure(merged, stats)
    write_report(stats)
    return stats


def make_figure(merged, stats):
    import matplotlib.pyplot as plt

    sf = {}
    for tag in SAMPLES:
        with open(os.path.join(S65, f"{tag}_scalefactors.json")) as f:
            sf[tag] = json.load(f)["tissue_lowres_scalef"]

    col_ldam = f"q50_{LDAM_LABEL}"
    vmax_a = max(merged[t][col_ldam].quantile(0.99) for t in SAMPLES)
    vmax_s = max(merged[t]["score_LDAM"].quantile(0.99) for t in SAMPLES)

    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    for i, tag in enumerate(["CTRL", "HS"]):
        j = merged[tag]
        img = plt.imread(os.path.join(S65, f"{tag}_tissue_lowres.png"))
        x = j["pxl_col_in_fullres"] * sf[tag]
        y = j["pxl_row_in_fullres"] * sf[tag]
        gsm = "GSM8778056 A18-43" if tag == "CTRL" else "GSM8778057 276A MTLE-HS"
        for k, (vals, vmax, title) in enumerate([
                (j[col_ldam], vmax_a,
                 f"cell2location {LDAM_LABEL} abundance (q50)"),
                (j["score_LDAM"], vmax_s, "LDAM gene-set score (scanpy)")]):
            ax = axes[i, k]
            ax.imshow(img, alpha=0.35)
            sc_ = ax.scatter(x, y, c=vals, s=14, cmap="viridis",
                             vmin=0, vmax=vmax, linewidths=0)
            ax.set_title(f"{tag} ({gsm})\n{title}", fontsize=9)
            ax.set_aspect("equal"); ax.invert_yaxis()
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            fig.colorbar(sc_, ax=ax, fraction=0.046, pad=0.02)
    mi = stats["morans_I"]
    fig.suptitle(
        "cell2location spatial mapping of LDAM-high microglia substate "
        f"(Moran's I: CTRL {mi['CTRL']['I']:.3f} (p={mi['CTRL']['perm_p']:.3f}), "
        f"HS {mi['HS']['I']:.3f} (p={mi['HS']['perm_p']:.3f}))",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fp = os.path.join(FIGDIR, "Fig_cell2location_maps.png")
    fig.savefig(fp, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log("figure saved " + fp)


def write_report(stats):
    refmeta = load_json(os.path.join(OUT, "reference_meta.json"))
    mapmeta = {t: load_json(os.path.join(OUT, f"mapping_meta_{t}.json"))
               for t in SAMPLES}
    hv = stats["HS_vs_CTRL"][f"q50_{LDAM_LABEL}"]
    cor = stats["correlation_with_geneset_scores"]
    lines = []
    A = lines.append
    A("# cell2location spatial deconvolution of LDAM-high microglia in MTLE-HS hippocampus")
    A("")
    A("## Aim")
    A("Quantify the spatial distribution of the LDAM (lipid-droplet-accumulating "
      "microglia) substate in human hippocampus Visium sections "
      "(GSE288888: GSM8778056 control, GSM8778057 MTLE-HS) by cell2location "
      "deconvolution, and test whether LDAM gene-set score hotspots coincide "
      "with mapped LDAM-high microglia.")
    A("")
    A("## Reference")
    A(f"- Source: {refmeta['ref_study']}, file `{os.path.basename(REF_H5AD)}`.")
    A("- Rationale: it is the only project reference with human temporal-lobe "
      "cells carrying both cell-type annotation and microglia sub-cluster "
      "structure with per-cell LDAM scores; no human hippocampus reference with "
      "LDAM substates exists in the workspace.")
    A(f"- Cells used: {refmeta['ref_cells_used']} (stratified subsample, "
      f"max {refmeta['max_cells_per_label']} per label); genes used: "
      f"{refmeta['ref_genes_used']} (cell2location filter_genes + MT-gene "
      "removal + intersection with Visium panel).")
    A("- Microglia were split into their leiden sub-clusters; sub-cluster 35 "
      "has the highest mean per-cell LDAM_score and was mapped as the "
      "`MG_LDAMhi` cell state; other microglia sub-clusters were kept as "
      "separate labels (`MG_c3/12/13/34`).")
    A("- Substate LDAM_score summary (reference):")
    for r in refmeta["microglia_substate_LDAM_score"]:
        A(f"  - {r['c2l_label']}: n={int(r['count'])}, "
          f"mean={r['mean']:.4f}, median={r['median']:.4f}")
    A(f"- RegressionModel (NB regression, batch=sample): max_epochs="
      f"{refmeta['ref_max_epochs']} (chunked, checkpointed every 50); "
      f"ELBO {refmeta['ref_elbo_first']:.4g} -> {refmeta['ref_elbo_last']:.4g}.")
    A("")
    A("## Mapping")
    for tag in ["CTRL", "HS"]:
        m = mapmeta[tag]
        A(f"- {tag}: {m['n_spots']} spots, {m['n_signature_genes']} signature "
          f"genes, N_cells_per_location={m['N_cells_per_location']}, "
          f"detection_alpha={m['detection_alpha']}, epochs={m['map_epochs']}, "
          f"ELBO {m['map_elbo_first']:.4g} -> {m['map_elbo_last']:.4g}.")
    A("")
    A("## Results")
    A("### LDAM-high microglia abundance, HS vs control")
    A(f"- median q50 abundance: HS {hv['median_HS']:.5f} vs CTRL "
      f"{hv['median_CTRL']:.5f}; mean {hv['mean_HS']:.5f} vs {hv['mean_CTRL']:.5f}.")
    A(f"- Mann-Whitney U={hv['mannwhitney_U']:.0f}, p={hv['p_value']:.3g}, "
      f"rank-biserial r={hv['rank_biserial_r_positive_means_HS_higher']:.3f} "
      "(positive = higher in HS).")
    tm = stats["HS_vs_CTRL"]["q50_total_microglia"]
    A(f"- Total microglia abundance: median HS {tm['median_HS']:.4f} vs CTRL "
      f"{tm['median_CTRL']:.4f}, rank-biserial r="
      f"{tm['rank_biserial_r_positive_means_HS_higher']:.3f}.")
    rv = stats["HS_vs_CTRL"][f"rel_{LDAM_LABEL}"]
    A(f"- Relative (compositional) abundance of {LDAM_LABEL} among mapped "
      f"cells: median HS {rv['median_HS']:.4f} vs CTRL {rv['median_CTRL']:.4f}, "
      f"Mann-Whitney p={rv['p_value']:.3g}, rank-biserial r="
      f"{rv['rank_biserial_r_positive_means_HS_higher']:.3f}.")
    A("- Note: absolute abundance is LOWER in HS than control while the LDAM "
      "gene-set score is HIGHER in HS (results/65); see Limitations for "
      "interpretation of this discrepancy.")
    A("")
    A("### Correlation with existing spot-level gene-set scores (results/65)")
    for k, d in cor.items():
        for tag, v in d.items():
            A(f"- {k} [{tag}]: Spearman rho={v['spearman_rho']:.3f}, "
              f"p={v['p_value']:.3g}, n={v['n']}")
    A("")
    A("### Spatial autocorrelation (Moran's I, array-grid neighbours)")
    for tag, v in stats["morans_I"].items():
        A(f"- {tag}: I={v['I']:.3f}, permutation p={v['perm_p']:.4f} "
          f"({v['n_perm']} permutations).")
    A("")
    A("### Hotspots (top quartile of LDAM-high microglia abundance)")
    for tag, v in stats["hotspots_top_quartile"].items():
        A(f"- {tag}: threshold(q75)={v['threshold_q75']:.5f}, "
          f"n={v['n_hotspots']} spots; median LDAM gene-set score in hotspots "
          f"{v['hotspot_median_score_LDAM']:.4f} vs "
          f"{v['nonhotspot_median_score_LDAM']:.4f} outside; region composition: "
          + ", ".join(f"{k} {vv:.2f}" for k, vv in
                      v["region_fraction_of_hotspots"].items()))
    A("")
    A("## Limitations")
    A("- Cross-slice direction discrepancy: absolute MG_LDAMhi abundance (and "
      "total mapped microglia) is lower in the HS section than in control, "
      "whereas the LDAM gene-set score is higher in HS. Possible causes: the "
      "reference is temporal cortex rather than hippocampus, the MG_LDAMhi "
      "label is an operational leiden sub-cluster, cell-density/detection "
      "differences between the two sections, and the sclerotic, neuron-depleted "
      "composition of MTLE-HS tissue. Within-section spatial agreement "
      "(CTRL Spearman rho=0.57 between MG_LDAMhi abundance and the LDAM "
      "gene-set score) is more trustworthy than cross-section absolute levels; "
      "compositional (relative) abundance is reported alongside.")
    A("- One Visium section per group (n=1 vs n=1): all HS-vs-control "
      "statistics are descriptive proof-of-concept, without biological "
      "replication; spot-level tests treat spots as units and do not estimate "
      "patient-level variance.")
    A("- The reference is human temporal-lobe (cortical) snRNA-seq from TLE "
      "surgical tissue, not hippocampus; hippocampus-specific populations "
      "(e.g. dentate granule cells) are absent from the reference, so their "
      "signal is absorbed by related types and absolute abundances are biased.")
    A("- The reference contains no curator-defined LDAM label: the `MG_LDAMhi` "
      "state is the leiden sub-cluster 35 with the highest per-cell LDAM_score "
      "in this dataset, i.e. an operational definition.")
    A("- CPU-only training: reference max_epochs capped at 250 and mapping "
      "epochs set by a wall-clock budget (~25 min/slice); ELBO trajectories "
      "above indicate convergence status.")
    A("- Rat spatial data were not deconvolved: the human reference does not "
      "match rat biology and cross-species cell2location mapping is unreliable.")
    A("")
    A("## Outputs")
    A("- `reference_signatures.csv` - per-label expression signatures "
      "(means_per_cluster_mu_fg).")
    A("- `abundance_per_spot.csv.gz` - q05/q50/q95 cell-state abundance per spot.")
    A("- `stats_summary.json` - all statistics reported here.")
    A("- `figures/v25/Fig_cell2location_maps.png` - spatial maps vs LDAM gene-set score.")
    with open(os.path.join(OUT, "c2l_report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    log("report written")


def main():
    inf_aver = run_reference()
    abund_all = {}
    for tag in ["CTRL", "HS"]:
        a, _ = map_sample(tag, inf_aver)
        abund_all[tag] = a
    pd.concat([clean_abund_columns(a) for a in abund_all.values()],
              ignore_index=True).to_csv(
        os.path.join(OUT, "abundance_per_spot.csv.gz"),
        index=False, compression="gzip")
    log("combined abundance table saved")
    downstream(abund_all)
    log("ALL DONE")


if __name__ == "__main__":
    main()
