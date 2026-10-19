"""Phase 2e: rewrite results/42_MR_coloc/summary.txt in UTF-8 with the full table,
merge cortex results into hub_genes_results.json."""
import json, os
import numpy as np

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
RES = os.path.join(BASE, "results", "42_MR_coloc")

lg_all = json.load(open(os.path.join(RES, "LGALS3_results.json")))          # hipp + cortex vs ALL
hub_hipp = json.load(open(os.path.join(RES, "hub_genes_results.json")))     # hipp vs ALL
hub_ctx = json.load(open(os.path.join(RES, "hub_genes_cortex_results.json")))  # cortex vs ALL
lg_foc = json.load(open(os.path.join(RES, "LGALS3_focal_results.json")))    # vs FOCAL

def or_ci(b, se):
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)

def mr_str(r):
    mr = r.get("MR:IVW (multiplicative RE)") or r.get("MR:Wald ratio")
    if not mr:
        return "insufficient instruments"
    name = "IVW" if "MR:IVW (multiplicative RE)" in r else "Wald"
    OR, lo, hi = or_ci(mr["b"], mr["se"])
    s = f"{name} OR={OR:.3f} ({lo:.3f}-{hi:.3f}), P={mr['p']:.3g}"
    eg = r.get("MR:MR-Egger")
    if eg:
        s += f"; Egger OR={np.exp(eg['b']):.3f}, P={eg['p']:.3g}, int-P={eg['intercept_p']:.3g}"
    return s

hub_hipp = [r for r in hub_hipp if "hippocampus" in r["tissue"]]
merged = hub_hipp + hub_ctx
with open(os.path.join(RES, "hub_genes_results.json"), "w", encoding="utf-8") as f:
    json.dump(merged, f, indent=1, default=float)

L = []
L.append("MR + colocalization: LPCAT2 / S100A10 / CD9 / SYNGR2 / LGALS3")
L.append("exposure = brain eQTL (MetaBrain), outcome = ILAE 2018 epilepsy GWAS")
L.append("=" * 78)
L.append("")
L.append("DATA SOURCES")
L.append("  Exposure: MetaBrain 2021-07-23 release (de Klein et al. 2023):")
L.append("    cortex-EUR-80PCs (N=2693), hippocampus-EUR-30PCs (N=168).")
L.append("    https://download.metabrain.nl/2021-07-23-release/2021-07-23-<tissue>-chr<N>.txt.gz")
L.append("    Downloaded this session: cortex chr1/chr12/chr16/chr17 (full files deleted")
L.append("    after extraction of per-gene cis rows; raw per-gene tsv kept).")
L.append("    LGALS3 cortex: local SMR_analysis/LGALS3_cortex_eqtl.txt (6,709 cis SNPs).")
L.append("  Outcomes: ILAE 2018 GWAS (ILAE Consortium on Complex Epilepsies, Nat Commun")
L.append("    2018, PMID 30531953), EBI GWAS Catalog FTP:")
L.append("    all epilepsy   GCST007343 all_epilepsy_METAL.gz   (15,212 cases / 29,677 controls)")
L.append("    focal epilepsy GCST007352 focal_epilepsy_METAL.gz (9,671 cases / 29,677 controls)")
L.append("    hg19 METAL; beta/se from Zscore, Freq1, Weight (per-SNP effective N).")
L.append("  Failed sources (documented): OpenGWAS JWT expired; GTEx v8 API returned empty/")
L.append("    connection errors; official MAGMA binary host unreachable (see 41_MAGMA summary).")
L.append("  Harmonization: rsID + effect-allele match; palindromic SNPs with MAF>0.42 dropped.")
L.append("  Instruments: eQTL P<5e-8, clumped r2<0.1, +/-250 kb (1000G EUR, N=503); relaxed to")
L.append("    P<1e-5 only when <3 instruments (flagged RELAXED). MR: IVW multiplicative RE,")
L.append("    Wald ratio (1 SNP), MR-Egger (>=3 SNPs). Coloc: coloc.abf re-implementation")
L.append("    (Giambartolomei 2014; p1=p2=1e-4, p12=1e-5, W=0.15/0.20), unit-tested; coloc")
L.append("    region = +/-500 kb around gene TSS (hg38), SNPs harmonized in both datasets.")
L.append("")
L.append("RESULTS TABLE (outcome: all epilepsy, GCST007343, unless noted)")
L.append(f"{'gene':8} {'tissue (N)':26} {'nHarm':>5} {'nIV':>3}  {'MR (OR per unit expression)':52} {'PP.H3':>6} {'PP.H4':>6}")
for r in lg_all + merged + lg_foc:
    tis = r["tissue"].replace("MetaBrain ", "").replace(" EUR", "")
    if "focal" in r.get("outcome", "all").lower() or r.get("outcome", "").startswith("ILAE 2018 focal"):
        tis += " [FOCAL]"
    co = r.get("coloc", {})
    L.append(f"{r['gene']:8} {tis:26.26} {r['n_harmonized']:>5} {r['n_instruments']:>3}  {mr_str(r):52.52} "
             f"{co.get('PP_H3', float('nan')):6.4f} {co.get('PP_H4', float('nan')):6.4f}")
L.append("SYNGR2   hippocampus (168)              0   0  not present in MetaBrain hippocampus release        -      -")
L.append("")
L.append("Instrument modes:")
for r in hub_ctx:
    L.append(f"  {r['gene']} cortex: {r['instrument_mode']}; IVs: {', '.join(r.get('instruments', []))}")
L.append("")
L.append("INTERPRETATION (honest)")
L.append("  * No hub gene shows a significant MR association in brain eQTL:")
L.append("    LPCAT2 cortex IVW P=0.772 (nIV=4); S100A10 P=0.793 (nIV=2); CD9 Wald P=0.553")
L.append("    (nIV=1); SYNGR2 P=0.849 (nIV=2, relaxed threshold).")
L.append("  * LGALS3 cortex vs all epilepsy: IVW OR=0.948 (0.902-0.996), P=0.0338 (nIV=2) -")
L.append("    nominally significant PROTECTIVE, but coloc PP.H4=0.0088 (PP.H1=0.949) shows the")
L.append("    eQTL signal does not colocalize with the GWAS signal; vs focal epilepsy:")
L.append("    OR=0.955 (0.901-1.011), P=0.115, PP.H4=0.0084.")
L.append("  * Every coloc analysis across the 5 genes gives PP.H4 < 0.01: no shared causal")
L.append("    variant between brain eQTL of these genes and epilepsy risk. PP.H0 dominates in")
L.append("    the underpowered hippocampus (N=168).")
L.append("  * Combined with prior deCODE pQTL results (results/15_MR_LGALS3_EPI_*.csv; coloc")
L.append("    PP.H4=0.024), single-gene genetic causality is unsupported; the genetic-level")
L.append("    evidence rests on the gene-set enrichment (results/41_MAGMA).")
L.append("")
L.append("FILES: exposure_*_{hippocampus,cortex}.csv(+_raw.tsv), harmonized_*.csv,")
L.append("  LGALS3_results.json, LGALS3_focal_results.json, hub_genes_results.json,")
L.append("  hub_genes_cortex_results.json")
L.append("Scripts: scripts/mr_coloc_utils.py, p2_lgals3.py, p2_hubs.py, p2c_cortex.py,")
L.append("  p2d_lgals3_focal.py, p2e_update_summary.py")
with open(os.path.join(RES, "summary.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L))
print("summary.txt (UTF-8) + merged hub_genes_results.json written")
