# -*- coding: utf-8 -*-
"""
p35_commot_spatial_lr.py
Spatially constrained ligand-receptor analysis with COMMOT (optimal-transport,
Cang et al., Nat Methods 2023) on Visium data:

  - Human hippocampus (GSE288888): MTLE-HS vs control (n=1 vs n=1).
  - Rat epileptogenesis time series (GSE301260): D5/D10/D20/D40 x SE/CTRL (n=2/group).

Core questions:
  a. Spatial distribution / domain enrichment of LGALS3 (galectin-3) signaling.
  b. Does the spatial extent ("range") and strength of LGALS3 signaling expand
     along epileptogenesis (D5 -> D40), SE vs CTRL?
  c. Spot-level coupling (Spearman) of LR signaling with LDAM / ferroptosis scores.
  d. Directionality: astrocyte-enriched domain -> microglia-enriched domain flow.

NOTE on the LGALS3 pairs: CellChatDB (shipped with commot 0.0.3) contains only
LGALS9 pairs, NOT LGALS3. LGALS3-receptor pairs supported by the literature
(LRP1, ITGB1, ITGA3, PTPRC/CD45, CD44, SLC3A2/CD98hc) are therefore appended as
custom database rows labelled pathway "GALECTIN3_LIT". The full CellChatDB
(secreted + contact + ECM) is screened alongside; the narrative focuses on Gal-3.

Statistical discipline: n=1/1 (human), n=2/group (rat). All spot-level tests are
descriptive; effect sizes are reported, p-values only where meaningful.

Run:  gsmap_env\\Scripts\\python.exe scripts\\p35_commot_spatial_lr.py \
          --stage human|rat|summarize|figs|all [--sample ID] [--force]
Outputs: results/70_commot_spatial/ , figures/v25/Fig_commot_*.png
"""
import argparse
import gzip
import json
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
from scipy import sparse, stats

import commot as cm

ROOT = Path(__file__).resolve().parents[1]
DATA_H = ROOT / "data" / "spatial" / "GSE288888"
DATA_R = ROOT / "data" / "spatial" / "GSE301260"
OUT = ROOT / "results" / "70_commot_spatial"
PER_SAMPLE = OUT / "per_sample"
FIGD = ROOT / "figures" / "v25"
for d in (OUT, PER_SAMPLE, FIGD):
    d.mkdir(parents=True, exist_ok=True)

HUMAN_SAMPLES = {
    "CTRL": "GSM8778056_A18-43",
    "HS": "GSM8778057_276A",
}

# Literature-supported galectin-3 (LGALS3) receptors — CellChatDB does NOT ship
# LGALS3 pairs (only LGALS9), so these are appended as custom rows.
GAL3_RECEPTORS = ["LRP1", "ITGB1", "ITGA3", "PTPRC", "CD44", "SLC3A2"]

# Core ligand set for the fast "core" mode: Gal-3 pairs + ferroptosis/lipid/
# myeloid-relevant signaling (SPP1, APP, PSAP, TGFb, TAM, CSF1/IL34, CX3CL1,
# TNF). ~26 CellChatDB pairs + 6 custom Gal-3 pairs.
CORE_LIGANDS = ["SPP1", "APP", "PSAP", "TGFB1", "GAS6", "PROS1", "CSF1",
                "IL34", "CX3CL1", "TNF"]

# COMMOT distance threshold: dis_thr = DIS_MULT * median nearest-neighbour
# distance. Visium spot pitch is 100 um centre-to-centre, so 3x ~= 300 um
# (three spot rings). Retries widen the radius if OT fails.
DIS_MULT = 3.0
DIS_MULT_RETRY = [3.0, 4.5, 6.0]
SPOT_DIAMETER_UM = 55.0  # Visium spot diameter -> pixel-to-um conversion

MICROGLIA_MARKERS = ["AIF1", "CX3CR1", "P2RY12", "TMEM119"]
ASTRO_MARKERS = ["GFAP", "AQP4"]


# ---------------------------------------------------------------- loading
def load_human(sample_id):
    pre = DATA_H / sample_id
    with gzip.open(f"{pre}_matrix.mtx.gz", "rb") as f:
        mat = scipy.io.mmread(f)
    X = sparse.csr_matrix(mat).T.tocsr()
    feats = pd.read_csv(f"{pre}_features.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol", "ftype"])
    bc = pd.read_csv(f"{pre}_barcodes.tsv.gz", header=None, names=["barcode"])
    pos = pd.read_csv(f"{pre}_tissue_positions.csv.gz")
    ad = sc.AnnData(X=X, obs=bc.set_index("barcode"), var=feats.set_index("symbol"))
    ad.var_names_make_unique()
    pos = pos.set_index("barcode")
    ad.obs = ad.obs.join(pos[["in_tissue", "array_row", "array_col",
                              "pxl_row_in_fullres", "pxl_col_in_fullres"]])
    ad = ad[ad.obs["in_tissue"] == 1].copy()
    ad.obsm["spatial"] = ad.obs[["pxl_col_in_fullres",
                                 "pxl_row_in_fullres"]].to_numpy(float)
    with gzip.open(f"{pre}_scalefactors_json.json.gz", "rt") as f:
        sf = json.load(f)
    return ad, sf


def load_rat(sample_id):
    pre = DATA_R / sample_id
    with gzip.open(f"{pre}_matrix.mtx.gz", "rb") as f:
        mat = scipy.io.mmread(f)
    X = sparse.csr_matrix(mat).T.tocsr()
    feats = pd.read_csv(f"{pre}_features.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol", "ftype"])
    bc = pd.read_csv(f"{pre}_barcodes.tsv.gz", header=None, names=["barcode"])
    pos = pd.read_csv(f"{pre}_tissue_positions_list.csv.gz", header=None,
                      names=["barcode", "in_tissue", "array_row", "array_col",
                             "pxl_col_in_fullres", "pxl_row_in_fullres"])
    ad = sc.AnnData(X=X, obs=bc.set_index("barcode"), var=feats.set_index("symbol"))
    ad.var_names_make_unique()
    pos = pos.set_index("barcode")
    ad.obs = ad.obs.join(pos[["in_tissue", "array_row", "array_col",
                              "pxl_row_in_fullres", "pxl_col_in_fullres"]])
    ad = ad[ad.obs["in_tissue"] == 1].copy()
    ad.obsm["spatial"] = ad.obs[["pxl_col_in_fullres",
                                 "pxl_row_in_fullres"]].to_numpy(float)
    with gzip.open(f"{pre}_scalefactors_json.json.gz", "rt") as f:
        sf = json.load(f)
    return ad, sf


def qc_normalize(ad):
    ad.var["mt"] = ad.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    return ad


def median_nn_dist(coords):
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=2).fit(coords)
    d, _ = nn.kneighbors(coords)
    return float(np.median(d[:, 1]))


def title_case(g):
    return g[0].upper() + g[1:].lower()


# ---------------------------------------------------------------- core COMMOT
def build_ligrec(ad, species, log, lr_set="full"):
    """CellChatDB filtered by expression + custom Gal-3 rows.

    lr_set='full': whole CellChatDB (all signaling types), expression-filtered.
    lr_set='core': Gal-3 pairs + CORE_LIGANDS pairs only (fast mode)."""
    db_species = "human" if species == "human" else "mouse"
    df = cm.pp.ligand_receptor_database(database="CellChat", species=db_species,
                                        signaling_type=None)
    df = df.iloc[:, :3].copy()
    df.columns = ["ligand", "receptor", "pathway"]
    n_db = len(df)
    if lr_set == "core":
        core_ligs = [g if species == "human" else title_case(g)
                     for g in CORE_LIGANDS]
        dff = df[df["ligand"].str.split("_").str[0].isin(core_ligs)].copy()
        # keep only pairs whose genes are detectable at all
        keep = []
        for _, r in dff.iterrows():
            genes = r["ligand"].split("_") + r["receptor"].split("_")
            keep.append(all(g in ad.var_names for g in genes))
        dff = dff[keep]
        log["n_pairs_db"] = int(n_db)
        log["n_pairs_after_expr_filter"] = int(len(dff))
    else:
        dff = cm.pp.filter_lr_database(df, ad, heteromeric=True,
                                       heteromeric_delimiter="_",
                                       heteromeric_rule="min",
                                       filter_criteria="min_cell_pct",
                                       min_cell_pct=0.05)
        dff = dff.iloc[:, :3].copy()
        dff.columns = ["ligand", "receptor", "pathway"]
        log["n_pairs_db"] = int(n_db)
        log["n_pairs_after_expr_filter"] = int(len(dff))
    # custom LGALS3 pairs
    lig = "LGALS3" if species == "human" else "Lgals3"
    rows = []
    for rec in GAL3_RECEPTORS:
        r = rec if species == "human" else title_case(rec)
        if lig in ad.var_names and r in ad.var_names:
            rows.append((lig, r, "GALECTIN3_LIT"))
    custom = pd.DataFrame(rows, columns=["ligand", "receptor", "pathway"])
    log["gal3_pairs_input"] = [f"{a}-{b}" for a, b, _ in rows]
    out = pd.concat([dff, custom], ignore_index=True).drop_duplicates(
        subset=["ligand", "receptor"]).reset_index(drop=True)
    log["n_pairs_run"] = int(len(out))
    return out


def run_commot_one(sample, species, treatment, timepoint, force=False,
                   lr_set="full"):
    """Per-sample COMMOT run; writes per-spot CSV + meta JSON.
    lr_set='core' writes *_core_* files and never overwrites full outputs."""
    sfx = "_core" if lr_set == "core" else ""
    meta_path = PER_SAMPLE / f"{sample}{sfx}_meta.json"
    csv_path = PER_SAMPLE / f"{sample}{sfx}_lr_signal_per_spot.csv"
    top_path = PER_SAMPLE / f"{sample}{sfx}_top_pairs.csv"
    if meta_path.exists() and csv_path.exists() and not force:
        print(f"[skip] {sample} {lr_set} (outputs exist)")
        return
    if lr_set == "core" and (PER_SAMPLE / f"{sample}_meta.json").exists():
        try:
            m = json.loads((PER_SAMPLE / f"{sample}_meta.json").read_text())
            if m.get("status") == "ok":
                print(f"[skip] {sample} core (full output exists)")
                return
        except Exception:  # noqa: BLE001
            pass
    t0 = time.time()
    log = {"sample": sample, "species": species, "treatment": treatment,
           "timepoint": timepoint, "lr_set": lr_set,
           "dis_mult_tried": [], "status": "ok"}

    if species == "human":
        ad, sf = load_human(sample)
    else:
        ad, sf = load_rat(sample)
    ad = qc_normalize(ad)
    log["n_spots"] = int(ad.n_obs)
    log["n_genes"] = int(ad.n_vars)
    spot_diam = float(sf.get("spot_diameter_fullres", np.nan))
    um_per_px = SPOT_DIAMETER_UM / spot_diam if spot_diam == spot_diam else np.nan
    log["spot_diameter_fullres_px"] = spot_diam
    log["um_per_px"] = um_per_px

    coords = ad.obsm["spatial"]
    d_nn = median_nn_dist(coords)
    log["median_nn_dist_px"] = d_nn
    log["median_nn_dist_um"] = d_nn * um_per_px

    df_ligrec = build_ligrec(ad, species, log, lr_set=lr_set)

    # ---- OT solve with retry on failure / empty Gal-3 output ----
    adata = None
    err = None
    for mult in DIS_MULT_RETRY:
        dis_thr = mult * d_nn
        log["dis_mult_tried"].append(mult)
        try:
            adata = ad.copy()
            cm.tl.spatial_communication(
                adata, database_name="cellchat", df_ligrec=df_ligrec,
                heteromeric=True, heteromeric_rule="min",
                heteromeric_delimiter="_", dis_thr=dis_thr,
                pathway_sum=False)
            ssum = adata.obsm["commot-cellchat-sum-sender"]
            gal3_cols = [c for c in ssum.columns if c[2:].startswith(
                "LGALS3" if species == "human" else "Lgals3")]
            total = float(ssum[gal3_cols].to_numpy().sum()) if gal3_cols else 0.0
            if gal3_cols and total > 0:
                log["dis_thr_px"] = float(dis_thr)
                log["dis_thr_um"] = float(dis_thr * um_per_px)
                log["dis_mult_used"] = mult
                err = None
                break
            err = f"empty Gal-3 signaling at dis_mult={mult}"
            print(f"  [retry] {sample}: {err}")
        except Exception as e:  # noqa: BLE001
            err = f"{type(e).__name__}: {e}"
            print(f"  [retry] {sample} dis_mult={mult}: {err}")
            adata = None
    if adata is None or err is not None:
        log["status"] = "failed"
        log["error"] = err
        meta_path.write_text(json.dumps(log, indent=2))
        print(f"[FAIL] {sample}: {err}")
        return

    db = "cellchat"
    ssum = adata.obsm[f"commot-{db}-sum-sender"]
    rsum = adata.obsm[f"commot-{db}-sum-receiver"]
    lig_prefix = "LGALS3-" if species == "human" else "Lgals3-"
    gal3_s_cols = [c for c in ssum.columns if c.startswith("s-" + lig_prefix)]
    gal3_r_cols = [c for c in rsum.columns if c.startswith("r-" + lig_prefix)]
    log["gal3_pairs_run"] = [c[2:] for c in gal3_s_cols]

    # ---- top pairs by total transported mass (full-db screen output) ----
    info_keys = [k for k in adata.obsp.keys()
                 if k.startswith(f"commot-{db}-") and "-" in k[len(f"commot-{db}-"):]]
    totals = []
    for k in info_keys:
        tag = k[len(f"commot-{db}-"):]
        if tag in ("total-total",):
            continue
        M = adata.obsp[k]
        totals.append((tag, float(M.sum())))
    top = pd.DataFrame(totals, columns=["lr_pair", "total_mass"]) \
        .sort_values("total_mass", ascending=False)
    top.to_csv(top_path, index=False)
    log["top20_pairs"] = top.head(20).to_dict("records")
    log["n_pairs_with_mass"] = int((top["total_mass"] > 0).sum())

    # ---- per-spot table ----
    out = pd.DataFrame({"barcode": adata.obs_names})
    out["pxl_col_in_fullres"] = coords[:, 0]
    out["pxl_row_in_fullres"] = coords[:, 1]
    out["array_row"] = adata.obs["array_row"].to_numpy()
    out["array_col"] = adata.obs["array_col"].to_numpy()
    for c in gal3_s_cols:
        out["gal3_sender_" + c[2 + len(lig_prefix):]] = ssum[c].to_numpy()
    for c in gal3_r_cols:
        out["gal3_receiver_" + c[2 + len(lig_prefix):]] = rsum[c].to_numpy()
    out["gal3_sender_total"] = out[[c for c in out.columns
                                    if c.startswith("gal3_sender_")]].sum(axis=1)
    out["gal3_receiver_total"] = out[[c for c in out.columns
                                      if c.startswith("gal3_receiver_")]].sum(axis=1)

    # ---- direction vector fields for Gal-3 pairs (k=5) ----
    vf_s = np.zeros((adata.n_obs, 2))
    vf_r = np.zeros((adata.n_obs, 2))
    n_vf = 0
    for pair in log["gal3_pairs_run"]:
        lig, rec = pair.split("-")
        try:
            cm.tl.communication_direction(adata, database_name=db,
                                          lr_pair=(lig, rec), k=5)
            ks = f"commot_sender_vf-{db}-{lig}-{rec}"
            kr = f"commot_receiver_vf-{db}-{lig}-{rec}"
            if ks in adata.obsm.keys():
                vf_s += np.asarray(adata.obsm[ks])
                vf_r += np.asarray(adata.obsm[kr])
                n_vf += 1
        except Exception as e:  # noqa: BLE001
            print(f"  [warn] direction {pair}: {e}")
    out["vf_sender_x"], out["vf_sender_y"] = vf_s[:, 0], vf_s[:, 1]
    out["vf_receiver_x"], out["vf_receiver_y"] = vf_r[:, 0], vf_r[:, 1]
    log["n_direction_fields"] = n_vf

    # ---- cell-type enrichment scores computed from expression ----
    def mean_expr(genes):
        genes = [g for g in genes if g in adata.var_names]
        if not genes:
            return np.full(adata.n_obs, np.nan)
        X = adata[:, genes].X
        if sparse.issparse(X):
            X = X.toarray()
        return X.mean(axis=1)
    if species == "human":
        out["astro_score"] = mean_expr(ASTRO_MARKERS)
        out["microglia_score"] = mean_expr(MICROGLIA_MARKERS)
    else:
        out["astro_score"] = mean_expr([title_case(g) for g in ASTRO_MARKERS])
        out["microglia_score"] = mean_expr([title_case(g) for g in MICROGLIA_MARKERS])

    out.to_csv(csv_path, index=False)
    log["runtime_s"] = round(time.time() - t0, 1)
    meta_path.write_text(json.dumps(log, indent=2))
    print(f"[done] {sample}: {log['n_spots']} spots, "
          f"{log['n_pairs_run']} pairs, {log['runtime_s']} s, "
          f"dis_thr={log['dis_thr_um']:.0f} um")


# ---------------------------------------------------------------- annotation merge
def attach_annotations():
    """Merge per-spot COMMOT output with GraphST domains / region labels / axis
    scores; writes lr_signal_per_spot.parquet (all samples).
    Prefers full-database outputs; falls back to core-set outputs per sample."""
    dom = pd.read_parquet(ROOT / "results" / "67_graphst" /
                          "domains_per_sample.parquet")
    dom["barcode_base"] = dom["barcode"].str.replace(r"-(CTRL|HS)$", "",
                                                     regex=True)
    # collect metas: full preferred over core per sample
    metas = {}
    for meta_path in sorted(PER_SAMPLE.glob("*_meta.json")):
        meta = json.loads(meta_path.read_text())
        if meta["status"] != "ok":
            continue
        sample = meta["sample"]
        lr_set = meta.get("lr_set", "full")
        if sample not in metas or lr_set == "full":
            metas[sample] = (meta_path, meta)
    meta_index = {}
    frames = []
    for sample, (meta_path, meta) in sorted(metas.items()):
        sfx = "_core" if meta.get("lr_set") == "core" else ""
        df = pd.read_csv(PER_SAMPLE / f"{sample}{sfx}_lr_signal_per_spot.csv")
        df["sample"] = sample
        df["species"] = meta["species"]
        df["treatment"] = meta["treatment"]
        df["timepoint"] = meta["timepoint"]
        df["lr_set"] = meta.get("lr_set", "full")
        meta_index[sample] = {"meta_file": meta_path.name,
                              "lr_set": meta.get("lr_set", "full"),
                              "um_per_px": meta.get("um_per_px"),
                              "dis_thr_um": meta.get("dis_thr_um"),
                              "n_pairs_run": meta.get("n_pairs_run")}
        df["barcode_base"] = df["barcode"].str.replace(r"-(CTRL|HS)$", "",
                                                       regex=True)
        d = dom[dom["sample"] == (
            f"human_{meta['treatment']}" if meta["species"] == "human"
            else sample)]
        d = d[["barcode_base", "graphst_domain", "graphst_annotation",
               "region", "score_LDAM", "score_Ferroptosis", "score_Module5",
               "score_Microglia"]].drop_duplicates("barcode_base")
        df = df.merge(d, on="barcode_base", how="left")
        frames.append(df)
    all_df = pd.concat(frames, ignore_index=True)
    all_df.to_parquet(OUT / "lr_signal_per_spot.parquet", index=False)
    all_df.to_csv(OUT / "lr_signal_per_spot.csv", index=False)
    (OUT / "sample_meta_index.json").write_text(json.dumps(meta_index,
                                                           indent=2))
    print(f"[merge] {all_df.shape} -> lr_signal_per_spot.parquet "
          f"({len(frames)} samples)")
    return all_df


# ---------------------------------------------------------------- summarize
def spearman(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 20:
        return np.nan, np.nan
    r, p = stats.spearmanr(a[m], b[m])
    return float(r), float(p)


def summarize():
    df = attach_annotations()
    df = df.copy()
    df["recv_pos"] = df["gal3_receiver_total"] > 0
    df["send_pos"] = df["gal3_sender_total"] > 0

    # ---------- per-domain LR signal density ----------
    rows = []
    for (sample, ann), g in df.groupby(["sample", "graphst_annotation"],
                                       dropna=False):
        rows.append({
            "sample": sample, "species": g["species"].iloc[0],
            "treatment": g["treatment"].iloc[0],
            "timepoint": g["timepoint"].iloc[0],
            "domain_annotation": ann,
            "n_spots": len(g),
            "gal3_sender_mean": g["gal3_sender_total"].mean(),
            "gal3_receiver_mean": g["gal3_receiver_total"].mean(),
            "frac_send_pos": g["send_pos"].mean(),
            "frac_recv_pos": g["recv_pos"].mean(),
            "score_LDAM_mean": g["score_LDAM"].mean(),
            "score_Ferroptosis_mean": g["score_Ferroptosis"].mean(),
        })
    dom_sum = pd.DataFrame(rows)
    dom_sum.to_csv(OUT / "domain_lr_summary.csv", index=False)

    # ---------- per-sample timecourse + correlations ----------
    tc = {}
    corr_rows = []
    for sample, g in df.groupby("sample"):
        sp, tr, tp = (g["species"].iloc[0], g["treatment"].iloc[0],
                      g["timepoint"].iloc[0])
        recv = g["gal3_receiver_total"].to_numpy()
        nz = recv[recv > 0]
        q90 = float(np.quantile(nz, 0.9)) if len(nz) else 0.0
        entry = {
            "species": sp, "treatment": tr, "timepoint": tp,
            "n_spots": int(len(g)),
            "gal3_sender_mean": float(g["gal3_sender_total"].mean()),
            "gal3_receiver_mean": float(recv.mean()),
            "gal3_receiver_mean_pos": float(nz.mean()) if len(nz) else 0.0,
            "frac_send_pos": float(g["send_pos"].mean()),
            "frac_recv_pos": float(g["recv_pos"].mean()),
            "frac_recv_top10pct": float((recv >= q90).mean()) if q90 > 0 else 0.0,
        }
        tc[sample] = entry
        for score in ["score_LDAM", "score_Ferroptosis"]:
            for sig in ["gal3_sender_total", "gal3_receiver_total"]:
                r, p = spearman(g[sig].to_numpy(float),
                                g[score].to_numpy(float))
                corr_rows.append({"sample": sample, "species": sp,
                                  "treatment": tr, "timepoint": tp,
                                  "signal": sig, "score": score,
                                  "spearman_r": r, "p": p})
    (OUT / "timecourse_lr.json").write_text(json.dumps(tc, indent=2))
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(OUT / "signal_score_correlations.csv", index=False)

    # ---------- directionality: astro-domain -> microglia-domain flow ----------
    meta_index = json.loads((OUT / "sample_meta_index.json").read_text())
    dir_rows = []
    for sample, g in df.groupby("sample"):
        ann = g["graphst_annotation"].fillna("unknown")
        astro_dom = g.groupby(ann)["astro_score"].mean().idxmax()
        mic_dom = g.groupby(ann)["microglia_score"].mean().idxmax()
        cA = g.loc[ann == astro_dom, ["pxl_col_in_fullres",
                                      "pxl_row_in_fullres"]].mean().to_numpy()
        cM = g.loc[ann == mic_dom, ["pxl_col_in_fullres",
                                    "pxl_row_in_fullres"]].mean().to_numpy()
        dvec = cM - cA
        nrm = np.linalg.norm(dvec)
        if nrm == 0 or astro_dom == mic_dom:
            continue
        u = dvec / nrm
        vf = g[["vf_sender_x", "vf_sender_y"]].to_numpy()
        proj = vf @ u
        send_pos = g["gal3_sender_total"] > 0
        fwd = proj[(ann == astro_dom) & send_pos]
        rev = proj[(ann == mic_dom) & send_pos]
        dir_rows.append({
            "sample": sample, "species": g["species"].iloc[0],
            "treatment": g["treatment"].iloc[0],
            "timepoint": g["timepoint"].iloc[0],
            "astro_domain": astro_dom, "microglia_domain": mic_dom,
            "domain_dist_um": float(nrm * meta_index[sample]["um_per_px"]),
            "flow_astro_to_mic_mean_proj": float(np.mean(fwd)) if len(fwd) else np.nan,
            "flow_mic_to_astro_mean_proj": float(np.mean(rev)) if len(rev) else np.nan,
            "n_sender_spots_astro": int(len(fwd)),
            "n_sender_spots_mic": int(len(rev)),
        })
    dir_df = pd.DataFrame(dir_rows)
    dir_df.to_csv(OUT / "direction_flow_summary.csv", index=False)

    # ---------- human HS vs CTRL effect sizes ----------
    hum = df[df["species"] == "human"]
    human_es = {}
    for sig in ["gal3_sender_total", "gal3_receiver_total"]:
        x = hum.loc[hum["treatment"] == "HS", sig].to_numpy()
        y = hum.loc[hum["treatment"] == "CTRL", sig].to_numpy()
        if len(x) and len(y):
            u_stat, p = stats.mannwhitneyu(x, y, alternative="two-sided")
            rbc = 2 * u_stat / (len(x) * len(y)) - 1
            human_es[sig] = {"HS_median": float(np.median(x)),
                             "CTRL_median": float(np.median(y)),
                             "rank_biserial": float(rbc), "p": float(p)}
    # HS vs CTRL frac recv pos
    fr = hum.groupby("treatment")["recv_pos"].mean().to_dict()
    human_es["frac_recv_pos"] = {k: float(v) for k, v in fr.items()}

    # ---------- rat timecourse summary (range expansion?) ----------
    rat = pd.DataFrame(tc).T
    rat = rat[rat["species"] == "rat"]
    rat_tc = {}
    for (tr, tp), g in rat.groupby(["treatment", "timepoint"]):
        rat_tc[f"{tr}_{tp}"] = {
            "n_samples": int(len(g)),
            "frac_recv_pos_mean": float(g["frac_recv_pos"].mean()),
            "frac_recv_pos_values": [float(v) for v in g["frac_recv_pos"]],
            "receiver_mean_pos_mean": float(g["gal3_receiver_mean_pos"].mean()),
            "receiver_mean_pos_values": [float(v)
                                         for v in g["gal3_receiver_mean_pos"]],
            "sender_mean_mean": float(g["gal3_sender_mean"].mean()),
            "sender_mean_values": [float(v) for v in g["gal3_sender_mean"]],
        }

    stats_summary = {
        "analysis": "COMMOT spatial ligand-receptor (optimal transport, "
                    "Cang et al. Nat Methods 2023)",
        "commot_version_note": "commot 0.0.3 (no __version__ attr; via "
                               "importlib.metadata)",
        "database": "CellChatDB (commot built-in; human for human samples, "
                    "mouse for rat) + custom LGALS3 literature pairs "
                    "(GALECTIN3_LIT: LRP1, ITGB1, ITGA3, PTPRC, CD44, SLC3A2)",
        "distance_threshold": f"dis_thr = {DIS_MULT}x median nearest-neighbour "
                              "distance (~300 um for Visium 100 um pitch); "
                              "retries at 4.5x/6x on failure",
        "n_design": {"human": "n=1 HS vs n=1 CTRL (descriptive only)",
                     "rat": "n=2 per timepoint x treatment (effect sizes only)"},
        "human_HS_vs_CTRL": human_es,
        "rat_timecourse": rat_tc,
        "correlation_table": "signal_score_correlations.csv",
        "direction_table": "direction_flow_summary.csv",
    }
    (OUT / "stats_summary.json").write_text(json.dumps(stats_summary, indent=2))
    print("[summarize] domain_lr_summary.csv / timecourse_lr.json / "
          "stats_summary.json written")
    return df


# ---------------------------------------------------------------- figures
def fig_human(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    hum = df[df["species"] == "human"]
    vmax = np.quantile(hum["gal3_receiver_total"], 0.99)
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    for ax, tr, title in zip(axes[0], ["CTRL", "HS"],
                             ["Control hippocampus (GSM8778056)",
                              "MTLE-HS hippocampus (GSM8778057)"]):
        g = hum[hum["treatment"] == tr]
        scat = ax.scatter(g["pxl_col_in_fullres"], g["pxl_row_in_fullres"],
                          c=g["gal3_receiver_total"], cmap="inferno",
                          s=8, vmin=0, vmax=vmax)
        ax.set_title(f"{title}\nGal-3 received signal (COMMOT)", fontsize=11)
        ax.invert_yaxis(); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(scat, ax=ax, shrink=0.7, label="receiver signal")
    # domain enrichment
    ax = axes[1, 0]
    ds = (hum.groupby(["graphst_annotation", "treatment"])
          ["gal3_receiver_total"].mean().unstack())
    ds = ds.sort_values("HS", ascending=False)
    x = np.arange(len(ds))
    ax.bar(x - 0.2, ds["CTRL"], 0.4, label="CTRL", color="#4C9BD6")
    ax.bar(x + 0.2, ds["HS"], 0.4, label="HS", color="#D64545")
    ax.set_xticks(x); ax.set_xticklabels(ds.index, rotation=35, ha="right",
                                         fontsize=9)
    ax.set_ylabel("mean Gal-3 receiver signal per spot")
    ax.set_title("Gal-3 signaling enrichment by domain (GraphST)", fontsize=11)
    ax.legend()
    # per-pair totals
    ax = axes[1, 1]
    pair_cols = [c for c in hum.columns if c.startswith("gal3_receiver_")
                 and c != "gal3_receiver_total"]
    pair_cols = [c for c in pair_cols if hum[c].fillna(0).sum() > 0]
    pt = hum.groupby("treatment")[pair_cols].sum().T
    pt.columns = ["CTRL", "HS"]
    pt.index = [c.replace("gal3_receiver_", "") for c in pt.index]
    pt = pt.sort_values("HS", ascending=True)
    ax.barh(np.arange(len(pt)) - 0.2, pt["CTRL"], 0.4, label="CTRL",
            color="#4C9BD6")
    ax.barh(np.arange(len(pt)) + 0.2, pt["HS"], 0.4, label="HS",
            color="#D64545")
    ax.set_yticks(np.arange(len(pt))); ax.set_yticklabels(pt.index, fontsize=9)
    ax.set_xlabel("summed receiver signal (all spots)")
    ax.set_title("Gal-3 receptor usage", fontsize=11)
    ax.legend()
    fig.suptitle("COMMOT spatial Gal-3 (LGALS3) signaling — human MTLE-HS vs "
                 "control hippocampus", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGD / "Fig_commot_human.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[fig] Fig_commot_human.png")


def fig_rat(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rat = df[df["species"] == "rat"]
    tp_order = ["D5", "D10", "D20", "D40"]
    tp_x = {t: i for i, t in enumerate(tp_order)}
    g = rat.groupby(["sample", "treatment", "timepoint"]).agg(
        frac_recv=("recv_pos", "mean"),
        recv_pos_mean=("gal3_receiver_total",
                       lambda v: v[v > 0].mean() if (v > 0).any() else 0.0),
        send_mean=("gal3_sender_total", "mean"),
    ).reset_index()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    metrics = [("frac_recv", "fraction of spots receiving Gal-3 signal",
                "Spatial extent of Gal-3 signaling"),
               ("recv_pos_mean", "mean receiver signal (positive spots)",
                "Gal-3 signal strength"),
               ("send_mean", "mean sender signal per spot",
                "Gal-3 sender activity")]
    colors = {"SE": "#D64545", "CTRL": "#4C9BD6"}
    for ax, (mcol, ylab, title) in zip(axes, metrics):
        for tr in ["CTRL", "SE"]:
            gg = g[g["treatment"] == tr]
            xs = [tp_x[t] for t in gg["timepoint"]]
            ax.scatter(xs, gg[mcol], color=colors[tr], s=45, alpha=0.85,
                       zorder=3, label=tr)
            med = gg.groupby("timepoint")[mcol].median().reindex(tp_order)
            ax.plot(range(4), med, color=colors[tr], lw=2, alpha=0.8)
        ax.set_xticks(range(4)); ax.set_xticklabels(tp_order)
        ax.set_xlabel("days after SE induction")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.3)
    axes[0].legend(title="treatment")
    fig.suptitle("Rat epileptogenesis time series — COMMOT Gal-3 (Lgals3) "
                 "signaling: extent and strength (points = samples, "
                 "n = 2/group; lines = medians)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIGD / "Fig_commot_rat_timecourse.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    print("[fig] Fig_commot_rat_timecourse.png")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["human", "rat", "summarize", "figs", "all"])
    ap.add_argument("--sample", default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--lr-set", default="full", choices=["full", "core"],
                    help="core = Gal-3 pairs + ferroptosis/lipid core ligands "
                         "only (fast); writes *_core_* files")
    args = ap.parse_args()
    sfx = "_core" if args.lr_set == "core" else ""

    if args.stage in ("human", "all"):
        for tr, sid in HUMAN_SAMPLES.items():
            if args.sample and sid != args.sample:
                continue
            try:
                run_commot_one(sid, "human", tr, "NA", force=args.force,
                               lr_set=args.lr_set)
            except Exception:  # noqa: BLE001
                traceback.print_exc()
                (PER_SAMPLE / f"{sid}{sfx}_meta.json").write_text(json.dumps(
                    {"sample": sid, "species": "human", "treatment": tr,
                     "lr_set": args.lr_set, "status": "failed",
                     "error": traceback.format_exc()[-2000:]}, indent=2))
    if args.stage in ("rat", "all"):
        qc = pd.read_csv(ROOT / "results" / "66_spatial_rat" / "sample_qc.csv")
        for _, r in qc.iterrows():
            if args.sample and r["sample"] != args.sample:
                continue
            try:
                run_commot_one(r["sample"], "rat", r["treatment"],
                               r["timepoint"], force=args.force,
                               lr_set=args.lr_set)
            except Exception:  # noqa: BLE001
                traceback.print_exc()
                (PER_SAMPLE / f"{r['sample']}{sfx}_meta.json").write_text(
                    json.dumps(
                        {"sample": r["sample"], "species": "rat",
                         "treatment": r["treatment"],
                         "timepoint": r["timepoint"], "lr_set": args.lr_set,
                         "status": "failed",
                         "error": traceback.format_exc()[-2000:]}, indent=2))
    df = None
    if args.stage in ("summarize", "figs", "all"):
        df = summarize()
    if args.stage in ("figs", "all"):
        fig_human(df)
        fig_rat(df)


if __name__ == "__main__":
    main()
