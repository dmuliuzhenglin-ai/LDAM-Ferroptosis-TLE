# -*- coding: utf-8 -*-
"""
p30_spatial_human.py
人脑海马空间转录组分析（10x Visium）：MTLE-HS vs 正常对照（GSE288888）
- GSM8778056_A18-43: 健康人对照海马 (CTRL)
- GSM8778057_276A:   耐药性 MTLE 伴海马硬化 (HS)

框架：n=1 vs n=1 proof-of-concept spatial mapping。
统计：spot 级分布差异（Mann-Whitney + rank-biserial 效应量）、
      spot 级 Spearman 相关、KNN 图 Moran's I（自写，squidpy 不可用）。

运行环境：gsmap_env\\Scripts\\python.exe scripts\\p30_spatial_human.py
输出：results/65_spatial_human/ 下的 spot 分数 CSV ×2、区域注释 CSV、统计 JSON、
      以及解压后的 lowres 组织图（供绘图脚本使用）。
"""
import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
from scipy import sparse, stats
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "spatial" / "GSE288888"
OUT = ROOT / "results" / "65_spatial_human"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLES = {
    "CTRL": "GSM8778056_A18-43",
    "HS": "GSM8778057_276A",
}

MICROGLIA_MARKERS = ["AIF1", "CX3CR1", "P2RY12", "TMEM119"]
SINGLE_GENES = ["LGALS3", "SPP1", "APOE", "FTH1", "FTL", "CH25H", "ACSL1", "CD9"]
REGION_MARKERS = {  # 大致亚区水平，不过度声称
    "DG": ["PROX1"],
    "CA1": ["NETO1"],
    "CA3": ["GRIK4"],  # GRIK4 在 DG/CA3 均有表达，仅作参考
    "astro": ["GFAP", "AQP4"],
    "neuron": ["SNAP25", "RBFOX3"],
}


def read_geneset(path):
    return [g.strip() for g in Path(path).read_text().splitlines() if g.strip()]


def load_sample(sample_id):
    """手动读取 Space Ranger MTX 三件套（v3 features，symbol 在第 2 列）。"""
    pre = DATA / sample_id
    with gzip.open(f"{pre}_matrix.mtx.gz", "rb") as f:
        mat = scipy.io.mmread(f)  # genes x spots
    X = sparse.csr_matrix(mat).T.tocsr()  # spots x genes
    feats = pd.read_csv(f"{pre}_features.tsv.gz", sep="\t", header=None,
                        names=["ensg", "symbol", "ftype"])
    bc = pd.read_csv(f"{pre}_barcodes.tsv.gz", header=None, names=["barcode"])
    pos = pd.read_csv(f"{pre}_tissue_positions.csv.gz")
    ad = sc.AnnData(X=X, obs=bc.set_index("barcode"),
                    var=feats.set_index("symbol"))
    ad.var_names_make_unique()
    pos = pos.set_index("barcode")
    ad.obs = ad.obs.join(pos[["in_tissue", "array_row", "array_col",
                              "pxl_row_in_fullres", "pxl_col_in_fullres"]])
    ad = ad[ad.obs["in_tissue"] == 1].copy()
    ad.obsm["spatial"] = ad.obs[["pxl_col_in_fullres",
                                 "pxl_row_in_fullres"]].to_numpy(float)
    return ad


def qc_normalize(ad):
    ad.var["mt"] = ad.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    return ad


def knn_morans_i(x, coords, k=6):
    """基于 KNN 空间权重（行标准化）的 Moran's I，置换检验 999 次。"""
    n = len(x)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()
    W = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    W = W.maximum(W.T)  # 对称化
    W = sparse.csr_matrix(W / W.sum(axis=1))
    S0 = W.sum()
    z = x - x.mean()
    I = (n / S0) * (z @ W @ z) / (z @ z)
    rng = np.random.default_rng(42)
    perm_I = np.empty(999)
    for i in range(999):
        zp = rng.permutation(z)
        perm_I[i] = (n / S0) * (zp @ W @ zp) / (zp @ zp)
    p = (np.sum(np.abs(perm_I) >= abs(I)) + 1) / (len(perm_I) + 1)
    return float(I), float(p)


def rank_biserial(x, y):
    """Mann-Whitney U + rank-biserial 相关效应量。x=HS, y=CTRL。"""
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    rbc = 2 * u / (len(x) * len(y)) - 1  # 正值 = HS 高于 CTRL
    return float(u), float(p), float(rbc)


def main():
    # ---------- 读数 + QC ----------
    ads = {}
    for tag, sid in SAMPLES.items():
        ad = qc_normalize(load_sample(sid))
        ad.obs["sample"] = tag
        ads[tag] = ad
        print(f"[load] {tag} {sid}: {ad.n_obs} in-tissue spots, {ad.n_vars} genes")

    # 解压 lowres 组织图供绘图用
    for tag, sid in SAMPLES.items():
        src = DATA / f"{sid}_tissue_lowres_image.png.gz"
        dst = OUT / f"{tag}_tissue_lowres.png"
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        sf = json.load(gzip.open(DATA / f"{sid}_scalefactors_json.json.gz", "rt"))
        json.dump(sf, open(OUT / f"{tag}_scalefactors.json", "w"))

    # ---------- 基因集 ----------
    ldam = read_geneset(ROOT / "genesets" / "LDAM_geneset.txt")
    ferro = read_geneset(ROOT / "genesets" / "Ferroptosis_geneset.txt")
    mod = pd.read_csv(ROOT / "results" / "59_metacell_wgcna" / "metacell_modules.csv")
    mod5 = mod.loc[mod["module"] == 5, "gene"].tolist()
    genesets = {"LDAM": ldam, "Ferroptosis": ferro, "Module5": mod5}

    # ---------- 合并打分（共同参照） ----------
    # 选择：两样本合并后用同一 control gene pool 打 score_genes(ctrl_size=50)，
    # 保证分数参照一致、可直接比较。
    combined = sc.concat(ads, label="batch", join="inner", index_unique="-")
    print(f"[merge] combined: {combined.n_obs} spots, {combined.n_vars} shared genes")

    avail = set(combined.var_names)
    for name, gs in genesets.items():
        gs_use = [g for g in gs if g in avail]
        print(f"[score] {name}: {len(gs_use)}/{len(gs)} genes detected")
        sc.tl.score_genes(combined, gs_use, ctrl_size=50,
                          score_name=f"score_{name}", random_state=0)

    # 单基因表达、微胶质富集度（标记均值）
    mg_use = [g for g in MICROGLIA_MARKERS if g in avail]
    mg_mat = combined[:, mg_use].X
    if sparse.issparse(mg_mat):
        mg_mat = mg_mat.toarray()
    combined.obs["microglia_enrichment"] = mg_mat.mean(axis=1)

    for g in SINGLE_GENES + ["PROX1", "NETO1", "GRIK4", "GFAP", "AQP4",
                             "SNAP25", "RBFOX3", "AIF1"]:
        if g in avail:
            v = combined[:, g].X
            combined.obs[f"expr_{g}"] = (v.toarray().ravel()
                                         if sparse.issparse(v) else v.ravel())

    # ---------- 区域注释 ----------
    # leidenalg 不可用，降级为 KMeans(k=8) on 表达 PCA（每样本各自聚类）。
    region_tables = []
    for tag in SAMPLES:
        sub = combined[combined.obs["sample"] == tag]
        Xp = sub.X.toarray() if sparse.issparse(sub.X) else sub.X
        n_pc = min(20, Xp.shape[1] - 1, Xp.shape[0] - 1)
        pcs = PCA(n_components=n_pc, random_state=0).fit_transform(Xp)
        km = KMeans(n_clusters=8, n_init=10, random_state=0).fit(pcs)
        sub.obs["cluster"] = km.labels_.astype(str)

        # 聚类 -> 亚区启发式注释：比较各 marker 组的簇均值
        marker_sets = {
            "DG-like (PROX1+)": ["expr_PROX1"],
            "CA1-like (NETO1+)": ["expr_NETO1"],
            "astro-enriched": ["expr_GFAP", "expr_AQP4"],
            "microglia-enriched": ["microglia_enrichment"],
        }
        rows = []
        cl_annot = {}
        # 先算各簇 marker 均值，再跨簇 z-score 比较，避免高本底 marker（GFAP）恒赢
        raw = {}
        for cl in sorted(sub.obs["cluster"].unique()):
            m = sub.obs["cluster"] == cl
            vals = {}
            for lab, cols in marker_sets.items():
                cols = [c for c in cols if c in sub.obs.columns]
                vals[lab] = float(sub.obs.loc[m, cols].to_numpy().mean()) if cols else 0.0
            raw[cl] = vals
        labs = list(marker_sets.keys())
        zmat = pd.DataFrame(raw).T[labs]  # cluster x marker-set
        zmat = (zmat - zmat.mean()) / zmat.std().replace(0, np.nan)
        for cl in sorted(sub.obs["cluster"].unique()):
            m = sub.obs["cluster"] == cl
            zr = zmat.loc[cl]
            best = zr.idxmax()
            # 最高 z-score 未超过 1（无显著富集）则标 mixed/other
            annot = best if (pd.notna(zr.max()) and zr.max() >= 1.0) else "mixed/other"
            cl_annot[cl] = annot
            rows.append({"sample": tag, "cluster": cl,
                         "n_spots": int(m.sum()), "annotation": annot,
                         **{f"mean_{k}": raw[cl][k] for k in labs},
                         **{f"z_{k}": float(zr[k]) for k in labs}})
        region_tables.append(pd.DataFrame(rows))
        # 同一样本内不同簇可能得到同一标签，加簇号区分
        combined.obs.loc[sub.obs_names, "region"] = [
            f"{cl_annot[c]}" for c in sub.obs["cluster"]]
        combined.obs.loc[sub.obs_names, "expr_cluster"] = sub.obs["cluster"].values
        print(f"[region] {tag}: " +
              "; ".join(f"c{r['cluster']}={r['annotation']}({r['n_spots']})"
                        for r in rows))

    region_df = pd.concat(region_tables, ignore_index=True)
    region_df.to_csv(OUT / "region_annotation.csv", index=False)

    # ---------- 保存 spot 分数 ----------
    score_cols = ["score_LDAM", "score_Ferroptosis", "score_Module5",
                  "microglia_enrichment", "region", "expr_cluster"]
    gene_cols = [c for c in combined.obs.columns if c.startswith("expr_")]
    for tag in SAMPLES:
        sub = combined.obs[combined.obs["sample"] == tag]
        df = sub[["array_row", "array_col", "pxl_row_in_fullres",
                  "pxl_col_in_fullres", "n_genes_by_counts",
                  "total_counts", "pct_counts_mt"] + score_cols + gene_cols].copy()
        df.insert(0, "barcode", df.index)
        df.to_csv(OUT / f"spot_scores_{tag}.csv", index=False)
        print(f"[save] spot_scores_{tag}.csv: {len(df)} spots")

    # ---------- 核心统计 ----------
    obs = combined.obs
    hs = obs[obs["sample"] == "HS"]
    ct = obs[obs["sample"] == "CTRL"]
    stats_out = {"framework": ("n=1 vs n=1 proof-of-concept spatial mapping; "
                               "spot-level distribution differences only, "
                               "no biological replicates"),
                 "n_spots": {t: int((obs['sample'] == t).sum()) for t in SAMPLES},
                 "scoring": ("scanpy tl.score_genes ctrl_size=50 on merged object "
                             "(common control-gene pool / common reference)"),
                 "comparisons": {}, "microglia_top_quartile": {},
                 "spearman": {}, "morans_i": {}}

    for ax in ["score_LDAM", "score_Ferroptosis", "score_Module5"]:
        u, p, rbc = rank_biserial(hs[ax].to_numpy(), ct[ax].to_numpy())
        stats_out["comparisons"][ax] = {
            "median_HS": float(hs[ax].median()), "median_CTRL": float(ct[ax].median()),
            "mean_HS": float(hs[ax].mean()), "mean_CTRL": float(ct[ax].mean()),
            "mannwhitney_U": u, "p_value": p, "rank_biserial_r": rbc}

    # 微胶质富集 spot（各样本内 AIF1 相关富集度 top 四分位）
    for tag, sub in [("HS", hs), ("CTRL", ct)]:
        thr = sub["microglia_enrichment"].quantile(0.75)
        mg = sub[sub["microglia_enrichment"] >= thr]
        for ax in ["score_LDAM", "score_Ferroptosis", "score_Module5"]:
            key = f"{ax}__{tag}"
            stats_out["microglia_top_quartile"][key] = {
                "n_spots": int(len(mg)),
                "median": float(mg[ax].median()), "mean": float(mg[ax].mean())}
    mg_hs = hs[hs["microglia_enrichment"] >= hs["microglia_enrichment"].quantile(0.75)]
    mg_ct = ct[ct["microglia_enrichment"] >= ct["microglia_enrichment"].quantile(0.75)]
    for ax in ["score_LDAM", "score_Ferroptosis", "score_Module5"]:
        u, p, rbc = rank_biserial(mg_hs[ax].to_numpy(), mg_ct[ax].to_numpy())
        stats_out["microglia_top_quartile"][f"{ax}__HS_vs_CTRL"] = {
            "mannwhitney_U": u, "p_value": p, "rank_biserial_r": rbc}

    # 轴分数 vs 微胶质富集度 Spearman（每样本分开）
    for tag, sub in [("HS", hs), ("CTRL", ct)]:
        for ax in ["score_LDAM", "score_Ferroptosis", "score_Module5"]:
            rho, p = stats.spearmanr(sub[ax], sub["microglia_enrichment"])
            stats_out["spearman"][f"{ax}__{tag}"] = {
                "rho": float(rho), "p_value": float(p), "n": int(len(sub))}

    # ---------- Moran's I（KNN k=6，每样本） ----------
    for tag in SAMPLES:
        sub = obs[obs["sample"] == tag]
        coords = sub[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
        for feat in ["score_LDAM", "score_Ferroptosis", "score_Module5",
                     "microglia_enrichment", "expr_LGALS3", "expr_SPP1",
                     "expr_APOE", "expr_FTH1"]:
            if feat not in sub.columns:
                continue
            I, p = knn_morans_i(sub[feat].to_numpy(float), coords)
            stats_out["morans_i"][f"{feat}__{tag}"] = {"I": I, "p_perm": p}
        print(f"[moran] {tag} done")

    # 亚区级汇总（若注释可靠）
    reg_stats = {}
    for tag in SAMPLES:
        sub = obs[obs["sample"] == tag]
        t = sub.groupby("region")[["score_LDAM", "score_Ferroptosis",
                                   "score_Module5"]].agg(["mean", "median", "count"])
        reg_stats[tag] = {r: {"n": int(t.loc[r, ("score_LDAM", "count")]),
                              "LDAM_mean": float(t.loc[r, ("score_LDAM", "mean")]),
                              "Ferroptosis_mean": float(t.loc[r, ("score_Ferroptosis", "mean")]),
                              "Module5_mean": float(t.loc[r, ("score_Module5", "mean")])}
                          for r in t.index}
    stats_out["region_summary"] = reg_stats

    with open(OUT / "stats_summary.json", "w") as f:
        json.dump(stats_out, f, indent=2, ensure_ascii=False)
    print(f"[save] stats_summary.json")
    print("DONE")


if __name__ == "__main__":
    main()
