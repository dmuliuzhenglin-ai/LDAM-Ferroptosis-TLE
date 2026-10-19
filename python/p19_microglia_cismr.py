"""Task: microglia cell-type-specific cis-MR (Bryois 2022 snRNA-seq eQTL, n=192)
for LDAM+ferroptosis axis core genes vs ILAE 2023 focal-HS / all / focal epilepsy.

Exposure: data/bryois/Microglia.<chr>.gz (5 cols, no header:
  <gene>_<ENSG> rsid pos P beta; beta w.r.t. effect_allele in snp_pos.txt.gz).
se derived: z = sign(beta)*Phi^-1(1-P/2); se = beta/z. F ~ z^2.
Instrument: lead SNP (min P) Wald ratio. Tier:
  P<5e-8 genome-wide significant; 5e-8<=P<1e-5 suggestive; P>=1e-5 not instrumentable.
Harmonisation: rsid join; flip vs GWAS Allele1 (effect allele); AT/GC palindromes dropped.
FDR: BH across all tested gene x outcome combos.
Outputs: results/64_microglia_cismr/
  microglia_cismr_results.csv, microglia_cismr_summary.json, comparison_with_blood.csv
"""
import os, json, gzip
import numpy as np, pandas as pd
from scipy import stats as st

BASE = os.getcwd()
OUT = os.path.join(BASE, "results", "64_microglia_cismr")
os.makedirs(OUT, exist_ok=True)

GENES = {1: ["S100A10"], 2: ["SYNGR2"], 3: ["PROK2"], 4: ["SCD5"],
         7: ["IKZF1", "JAZF1"], 12: ["CD9", "ETV6"], 14: ["LGALS3"],
         16: ["IRF8", "LPCAT2"], 21: ["RUNX1"]}
ALL_GENES = [g for v in GENES.values() for g in v]

# ---------------- exposure: stream chromosome files, keep target genes ----------------
def extract_gene_rows():
    rows = []
    for chrom, genes in GENES.items():
        fn = os.path.join(BASE, "data", "bryois", f"Microglia.{chrom}.gz")
        want = set(genes)
        with gzip.open(fn, "rt") as f:
            for line in f:
                head = line.split(" ", 1)[0]
                sym = head.split("_")[0]
                if sym in want:
                    parts = line.split()
                    rows.append(dict(gene=sym, gid=head, rsid=parts[1],
                                     p_eqtl=float(parts[3]), beta_x=float(parts[4]),
                                     chrom=chrom))
        print(f"chr{chrom} scanned", flush=True)
    return pd.DataFrame(rows)

# ---------------- snp annotation (stream-filter to needed rsids) ----------------
def load_snp_pos(need):
    keep = []
    fn = os.path.join(BASE, "data", "bryois", "snp_pos.txt.gz")
    for chunk in pd.read_csv(fn, sep="\t", chunksize=500_000, dtype=str):
        sub = chunk[chunk.SNP.isin(need)]
        if len(sub):
            keep.append(sub)
    sp = pd.concat(keep, ignore_index=True)
    sp["effect_allele"] = sp["effect_allele"].str.upper()
    sp["other_allele"] = sp["other_allele"].str.upper()
    return sp[["SNP", "effect_allele", "other_allele", "MAF"]].rename(columns={"SNP": "rsid"})

# ---------------- outcome loader (verbatim from p18_axis_cismr.py) ----------------
def load_outcome(path):
    raw = pd.read_csv(path, sep="\t", header=None, skiprows=1, dtype=str)
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
    """exp: rsid, ea_x (eQTL effect allele), nea_x; out: rsid, ea_y=Allele1(effect), nea_y.
    Returns harmonised frame + n palindromic dropped."""
    m = exp.merge(out, on="rsid", suffixes=("_x", "_y"))
    pal = (m.ea_x + m.nea_x).isin(["AT", "TA", "CG", "GC"])
    n_pal = int(pal.sum())
    m = m[~pal]
    sign = np.where((m.ea_x == m.ea_y) & (m.nea_x == m.nea_y), 1,
            np.where((m.ea_x == m.nea_y) & (m.nea_x == m.ea_y), -1,
            np.where((flip(m.ea_x) == m.ea_y) & (flip(m.nea_x) == m.nea_y), 1,
            np.where((flip(m.ea_x) == m.nea_y) & (flip(m.nea_x) == m.ea_y), -1, 0))))
    m = m[sign != 0].copy()
    m["sign"] = sign[sign != 0]
    m["beta_y"] = m["beta_y"] * m["sign"]
    return m, n_pal

def main():
    eq = extract_gene_rows()
    print("exposure rows:", len(eq), "genes found:", sorted(eq.gene.unique()), flush=True)
    tested = sorted(eq.gene.unique())
    not_tested = [g for g in ALL_GENES if g not in tested]

    sp = load_snp_pos(set(eq.rsid))
    eq = eq.merge(sp, on="rsid", how="inner")
    print("rows with allele annotation:", len(eq), flush=True)

    # derive z, se, F
    z_abs = st.norm.isf(eq.p_eqtl / 2.0)
    eq["z_x"] = np.sign(eq.beta_x) * z_abs
    eq["se_x"] = eq.beta_x / eq.z_x
    eq["F"] = eq.z_x ** 2

    def tier(p):
        return ("genome-wide significant instrument" if p < 5e-8 else
                "suggestive instrument" if p < 1e-5 else "not instrumentable (lead P>=1e-5)")

    outcomes = {
        "ILAE3_focalHS": load_outcome(os.path.join(BASE, "data", "gwas", "ILAE2023", "ILAE3_focal_HS_final.tbl")),
        "ILAE3_all": load_outcome(os.path.join(BASE, "data", "gwas", "ILAE2023", "ILAE3_Caucasian_all_epilepsy_final.tbl")),
        "ILAE3_focal": load_outcome(os.path.join(BASE, "data", "gwas", "ILAE2023", "ILAE3_Caucasian_focal_epilepsy_final.tbl")),
    }

    # lead SNP per gene: min eQTL P among SNPs harmonisable with the PRIMARY outcome
    # (focal-HS). If the strict min-P SNP is unavailable/unusable in the outcome GWAS,
    # fall back to the next-best harmonisable SNP -- threshold policy unchanged.
    exp_all = eq.rename(columns={"effect_allele": "ea_x", "other_allele": "nea_x"})
    hm_primary, n_pal_primary = harmonise(
        exp_all, outcomes["ILAE3_focalHS"].rename(
            columns={"ea": "ea_y", "nea": "nea_y", "beta": "beta_y", "se": "se_y"}))
    strict_lead = eq.sort_values("p_eqtl").groupby("gene").head(1).set_index("gene")
    lead = hm_primary.sort_values("p_eqtl").groupby("gene").head(1).set_index("gene")
    # prefer the strict min-P SNP when it is itself harmonisable (LD ties in P)
    for g in lead.index:
        s = strict_lead.loc[g]
        if len(hm_primary[(hm_primary.gene == g) & (hm_primary.rsid == s.rsid)]) == 1:
            lead.loc[g] = hm_primary[(hm_primary.gene == g) & (hm_primary.rsid == s.rsid)].iloc[0]
    inst_info = {}
    for g in tested:
        r = lead.loc[g]
        s = strict_lead.loc[g]
        inst_info[g] = dict(
            rsid=r.rsid, ea=r.ea_x, nea=r.nea_x,
            beta_x=r.beta_x, se_x=r.se_x, p_eqtl=r.p_eqtl, F=r.F,
            tier=tier(r.p_eqtl), maf=float(r.MAF),
            n_cis_snps_tested=int((eq.gene == g).sum()),
            strict_lead_rsid=s.rsid, strict_lead_p=float(s.p_eqtl),
            lead_selection=("min P overall" if s.rsid == r.rsid else
                            "min P among harmonisable (strict lead unavailable/unusable in focal-HS GWAS)"))

    res, harm_log = [], {}
    for oname, ogw in outcomes.items():
        ogw = ogw.rename(columns={"ea": "ea_y", "nea": "nea_y", "beta": "beta_y", "se": "se_y"})
        hm, n_pal = harmonise(exp_all, ogw)
        harm_log[oname] = dict(palindromic_dropped=n_pal, harmonised_rows=len(hm))
        print(oname, "palindromic dropped:", n_pal, "harmonised rows:", len(hm), flush=True)
        for g in tested:
            info = inst_info[g]
            sub = hm[(hm.gene == g) & (hm.rsid == info["rsid"])]
            row = dict(gene=g, outcome=oname, lead_rsid=info["rsid"],
                       ea=info["ea"], nea=info["nea"], maf=info["maf"],
                       beta_eqtl=info["beta_x"], se_eqtl=info["se_x"], p_eqtl=info["p_eqtl"],
                       F=info["F"], instrument_tier=info["tier"],
                       lead_selection=info["lead_selection"],
                       n_cis_snps_tested=info["n_cis_snps_tested"])
            if len(sub) == 1:
                r = sub.iloc[0]
                b = r.beta_y / r.beta_x
                se = r.se_y / abs(r.beta_x)
                p = 2 * st.norm.sf(abs(b / se))
                row.update(beta_outcome_aligned=r.beta_y, se_outcome=r.se_y,
                           beta_mr=b, se_mr=se, p_mr=p, harmonised=True)
            else:
                row.update(beta_outcome_aligned=np.nan, se_outcome=np.nan,
                           beta_mr=np.nan, se_mr=np.nan, p_mr=np.nan, harmonised=False)
            res.append(row)
    R = pd.DataFrame(res)
    # BH FDR across all tested gene x outcome combos
    mask = R.p_mr.notna()
    R["fdr"] = np.nan
    R.loc[mask, "fdr"] = st.false_discovery_control(R.loc[mask, "p_mr"])
    R = R.sort_values(["gene", "outcome"])
    R.to_csv(os.path.join(OUT, "microglia_cismr_results.csv"), index=False)

    # ---- comparison with eQTLGen whole blood (results/60_axis_cismr) ----
    blood_path = os.path.join(BASE, "results", "60_axis_cismr", "axis_cismr_results.csv")
    blood = pd.read_csv(blood_path)
    blood = blood[blood.gene.isin(ALL_GENES)][
        ["outcome", "gene", "n_inst", "lead_rsid", "beta_wald", "se_wald", "p_wald", "fdr_wald"]]
    blood = blood.rename(columns={c: "blood_" + c for c in
                                  ["n_inst", "lead_rsid", "beta_wald", "se_wald", "p_wald", "fdr_wald"]})
    grid = pd.MultiIndex.from_product([ALL_GENES, list(outcomes)], names=["gene", "outcome"]).to_frame(index=False)
    comp = grid.merge(R[["gene", "outcome", "lead_rsid", "instrument_tier", "F",
                         "beta_mr", "se_mr", "p_mr", "fdr", "harmonised"]],
                      on=["gene", "outcome"], how="left")
    comp = comp.merge(blood, on=["gene", "outcome"], how="left")
    comp["microglia_tested"] = comp.gene.isin(tested)
    comp.to_csv(os.path.join(OUT, "comparison_with_blood.csv"), index=False)

    # ---- summary json ----
    summary = dict(
        exposure="Bryois 2022 microglia snRNA-seq cis-eQTL (n=192), Zenodo 10.5281/zenodo.7276971",
        outcomes=list(outcomes),
        genes_target=ALL_GENES, genes_tested=tested, genes_not_tested=not_tested,
        not_tested_reason="gene absent from Bryois microglia cis-eQTL files (not expressed/tested in microglia panel); chr14 file re-downloaded from Zenodo (local copy was truncated) and LGALS3 confirmed absent",
        instruments=inst_info,
        harmonisation=harm_log,
        n_tests=int(mask.sum()),
        n_fdr_significant=int((R.fdr < 0.05).sum()),
        n_nominal_significant=int((R.p_mr < 0.05).sum()),
        min_p_per_outcome=R.loc[mask].groupby("outcome").p_mr.min().to_dict(),
    )
    json.dump(summary, open(os.path.join(OUT, "microglia_cismr_summary.json"), "w"),
              indent=2, default=float)
    print(json.dumps({k: v for k, v in summary.items() if k != "instruments"}, indent=2, default=str))
    print("== instruments ==")
    for g, i in inst_info.items():
        print(f"{g}: {i['rsid']} P={i['p_eqtl']:.3g} beta={i['beta_x']:.4f} se={i['se_x']:.4f} F={i['F']:.1f} [{i['tier']}]")
    print("== top MR results ==")
    print(R.loc[mask].sort_values("p_mr").head(12)[
        ["gene", "outcome", "lead_rsid", "beta_mr", "se_mr", "p_mr", "fdr"]].to_string(index=False))

if __name__ == "__main__":
    main()
