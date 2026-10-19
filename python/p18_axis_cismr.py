"""Task 7: systematic cis-MR screen of LDAM + ferroptosis axis genes (eQTLGen blood)
against ILAE 2023 focal-HS and ILAE 2018 all-epilepsy GWAS.
Instruments: OpenGWAS tophits p<5e-8 (batch API); lead-SNP Wald primary,
up-to-3-SNP (>=1Mb apart) IVW sensitivity. Palindromic SNPs dropped.
Env: OPENGWAS_JWT must be set.
Outputs: results/60_axis_cismr/
"""
import os, json, time
import numpy as np, pandas as pd
import requests
from scipy import stats as st

BASE = os.getcwd()
OUT = os.path.join(BASE, "results", "60_axis_cismr")
os.makedirs(OUT, exist_ok=True)
TOK = os.environ["OPENGWAS_JWT"]
H = {"Authorization": f"Bearer {TOK}"}

def mygene_map(genes):
    res = {}
    for i in range(0, len(genes), 1000):
        chunk = genes[i:i + 1000]
        r = requests.post("https://mygene.info/v3/query",
                          data={"q": ",".join(chunk), "scopes": "symbol",
                                "fields": "ensembl.gene", "species": "human"}, timeout=120)
        for x in r.json():
            if x.get("notfound"):
                continue
            e = x.get("ensembl")
            if isinstance(e, list):
                e = e[0] if e else None
            if e and "gene" in e:
                res[x["query"]] = e["gene"]
        time.sleep(0.5)
    return res

def tophits_batch(ids):
    out = []
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        for attempt in range(4):
            try:
                r = requests.post("https://api.opengwas.io/api/tophits",
                                  json={"id": chunk, "pval": 5e-8, "clump": 0},
                                  headers=H, timeout=180)
                if r.status_code == 200:
                    out.extend(r.json())
                    break
                time.sleep(10 * (attempt + 1))
            except Exception:
                time.sleep(10 * (attempt + 1))
        time.sleep(1.2)
    return out

def load_outcome(path, tag):
    raw = pd.read_csv(path, sep="\t", header=None, skiprows=1, dtype=str)
    # ILAE3 tbl: 11 tab fields; the last one holds "Direction Beta SE" space-joined
    raw.columns = ["CHR", "BP", "rsid", "ea", "nea", "Freq1", "FreqSE",
                   "Effective_N", "Z", "P", "tail"]
    parts = raw["tail"].str.strip().str.split(r"\s+", expand=True)
    raw["beta"] = pd.to_numeric(parts[1], errors="coerce")
    raw["se"] = pd.to_numeric(parts[2], errors="coerce")
    raw["ea"] = raw["ea"].str.upper(); raw["nea"] = raw["nea"].str.upper()
    return raw[["rsid", "ea", "nea", "beta", "se"]].dropna()

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
def flip(a):
    return "".join(COMP[x] for x in a)

def harmonise(exp, out):
    m = exp.merge(out, on="rsid", suffixes=("_x", "_y"))
    m = m[~((m.ea_x + m.nea_x).isin(["AT", "TA", "CG", "GC"]))]
    sign = np.where((m.ea_x == m.ea_y) & (m.nea_x == m.nea_y), 1,
            np.where((m.ea_x == m.nea_y) & (m.nea_x == m.ea_y), -1,
            np.where((flip(m.ea_x) == m.ea_y) & (flip(m.nea_x) == m.nea_y), 1,
            np.where((flip(m.ea_x) == m.nea_y) & (flip(m.nea_x) == m.ea_y), -1, 0))))
    m = m[sign != 0].copy()
    m["beta_y"] = m["beta_y"] * sign[sign != 0]
    return m

def main():
    ldam = set(pd.read_csv(os.path.join(BASE, "genesets", "LDAM_geneset.txt"), header=None)[0].str.strip())
    fer = set(pd.read_csv(os.path.join(BASE, "genesets", "Ferroptosis_geneset.txt"), header=None)[0].str.strip())
    genes = sorted(ldam | fer)
    print("axis genes:", len(genes), flush=True)
    cache = os.path.join(OUT, "symbol2ensg.json")
    if os.path.exists(cache):
        s2e = json.load(open(cache))
    else:
        s2e = mygene_map(genes)
        json.dump(s2e, open(cache, "w"))
    info = json.load(open(os.path.join(BASE, "results", "54_hub_pQTL_MR", "gwasinfo_all.json")))
    rows = [{"gene": g, "ensg": s2e[g], "gid": f"eqtl-a-{s2e[g]}"}
            for g in genes if g in s2e and f"eqtl-a-{s2e[g]}" in info]
    gm = pd.DataFrame(rows)
    print("with eQTLGen id:", len(gm), flush=True)
    th_cache = os.path.join(OUT, "tophits_raw.json")
    if os.path.exists(th_cache):
        th = json.load(open(th_cache))
    else:
        th = tophits_batch(gm.gid.tolist())
        json.dump(th, open(th_cache, "w"))
    th = pd.DataFrame(th)
    print("instrument rows:", len(th), flush=True)
    th = th[th.p < 5e-8]
    # select instruments per gene: lead + up to 2 more >=1Mb apart
    inst = []
    for gid, grp in th.groupby("id"):
        grp = grp.sort_values("p")
        sel = []
        for _, r in grp.iterrows():
            if all(abs(int(r["position"]) - int(s["position"])) >= 1_000_000 for s in sel):
                sel.append(r)
            if len(sel) == 3:
                break
        for s in sel:
            inst.append({"gid": gid, "rsid": s["rsid"], "ea": s["ea"].upper(), "nea": s["nea"].upper(),
                         "beta_x": s["beta"], "se_x": s["se"], "pos": int(s["position"])})
    inst = pd.DataFrame(inst)
    inst = inst.merge(gm, on="gid")
    inst.to_csv(os.path.join(OUT, "instruments.csv"), index=False)
    print("genes with instruments:", inst.gene.nunique(), flush=True)

    outcomes = {
        "ILAE3_focalHS": load_outcome(os.path.join(BASE, "data", "gwas", "ILAE2023", "ILAE3_focal_HS_final.tbl"), "fhs"),
        "ILAE3_all": load_outcome(os.path.join(BASE, "data", "gwas", "ILAE2023", "ILAE3_Caucasian_all_epilepsy_final.tbl"), "all"),
        "ILAE3_focal": load_outcome(os.path.join(BASE, "data", "gwas", "ILAE2023", "ILAE3_Caucasian_focal_epilepsy_final.tbl"), "foc"),
    }
    res = []
    for oname, ogw in outcomes.items():
        hm = harmonise(inst.rename(columns={"ea": "ea_x", "nea": "nea_x"}), ogw.rename(columns={"ea": "ea_y", "nea": "nea_y", "beta": "beta_y", "se": "se_y"}))
        for gene, grp in hm.groupby("gene"):
            grp = grp.sort_values("pos")
            # lead SNP Wald
            lead = grp.loc[grp.beta_x.abs().idxmax()]
            b = lead.beta_y / lead.beta_x
            se = lead.se_y / abs(lead.beta_x)
            p = 2 * st.norm.sf(abs(b / se))
            # IVW if >1 instrument
            piv, seiv = np.nan, np.nan
            if len(grp) > 1:
                bi = grp.beta_y / grp.beta_x
                vi = (grp.se_y / grp.beta_x) ** 2
                biv = np.sum(bi / vi) / np.sum(1 / vi)
                seiv = np.sqrt(1 / np.sum(1 / vi))
                piv = 2 * st.norm.sf(abs(biv / seiv))
            res.append(dict(outcome=oname, gene=gene, n_inst=len(grp),
                            lead_rsid=lead.rsid, beta_wald=b, se_wald=se, p_wald=p,
                            beta_ivw=biv if len(grp) > 1 else np.nan,
                            se_ivw=seiv, p_ivw=piv))
    R = pd.DataFrame(res)
    R["fdr_wald"] = R.groupby("outcome").p_wald.transform(
        lambda p: pd.Series(st.false_discovery_control(p.fillna(1)), index=p.index) if hasattr(st, "false_discovery_control") else p)
    R.to_csv(os.path.join(OUT, "axis_cismr_results.csv"), index=False)
    for oname in outcomes:
        sub = R[R.outcome == oname].sort_values("p_wald")
        print("==", oname, "tested:", len(sub), "FDR<0.05:", (sub.fdr_wald < 0.05).sum())
        print(sub.head(10)[["gene", "n_inst", "beta_wald", "p_wald", "fdr_wald"]].to_string(index=False))

if __name__ == "__main__":
    main()
