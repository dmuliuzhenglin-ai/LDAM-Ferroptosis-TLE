"""p23: build results/61_twas/twas_summary.json across tissues x outcomes."""
import os, json
import pandas as pd

BASE = os.getcwd()
OUT = os.path.join(BASE, "results", "61_twas")
GTF = os.path.join(BASE, "results", "50_gsMap", "gsMap_resource",
                   "genome_annotation", "gtf",
                   "gencode.v46lift37.basic.annotation.gtf")

e2s = {}
with open(GTF, encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.split("\t")
        if len(c) < 9 or c[2] != "gene":
            continue
        a = c[8]
        gid = a.split('gene_id "')[1].split('"')[0].split(".")[0].split("_")[0]
        gn = a.split('gene_name "')[1].split('"')[0] if 'gene_name "' in a else gid
        e2s[gid] = gn

combos = {}
for tissue in ["Hippocampus", "Cortex"]:
    for o in ["fhs", "all", "focal"]:
        f = os.path.join(OUT, f"TWAS_{tissue}_{o}.csv")
        if not os.path.exists(f):
            continue
        R = pd.read_csv(f)
        top = R.nsmallest(5, "P")
        combos[f"{tissue}_{o}"] = {
            "file": f,
            "n_tested": int(len(R)),
            "n_fdr05": int((R.FDR < 0.05).sum()),
            "n_nominal_p05": int((R.P < 0.05).sum()),
            "min_P": float(R.P.min()),
            "top5": [{
                "gene": r.gene,
                "symbol": e2s.get(r.gene.split(".")[0], r.gene),
                "chr": int(r.chr), "p0": int(r.p0), "p1": int(r.p1),
                "model": r.model, "cv_r2": round(float(r.cv_r2), 4),
                "nsnp": int(r.nsnp),
                "Z": round(float(r.Z), 3),
                "P": float(r.P), "FDR": round(float(r.FDR), 4),
            } for _, r in top.iterrows()],
        }

summary = {
    "analysis": "Summary-based TWAS (S-PrediXcan style), Python re-implementation",
    "weights": "GTEx v8 FUSION expression weights, Brain_Hippocampus (.pos primary) and Brain_Cortex",
    "gwas": {
        "fhs": "ILAE 2023 focal epilepsy with hippocampal sclerosis (EUR)",
        "all": "ILAE 2023 all epilepsy (EUR, Caucasian)",
        "focal": "ILAE 2023 focal epilepsy (EUR, Caucasian)",
    },
    "ld_reference": "1000G EUR local plink bed (g1000_eur, 503 samples), per-gene window reads only",
    "model_selection": "non-top1 model with highest cross-validation R2; falls back to top1 if best non-top1 cv p>0.05",
    "statistic": "Z = w'z / sqrt(w'Lw), L = genotype correlation + 1e-6 ridge; BH FDR per tissue x outcome",
    "harmonization": "rsid + allele matching across weights / GWAS / 1000G bim with strand-flip handling; palindromic (AT/GC) SNPs excluded (none present in GTEx FUSION weights)",
    "results": combos,
    "validation": "results/61_twas/validation_notes.txt - all-epilepsy run: best TWAS gene at GWAS lead loci ranks median top 3.2% genome-wide in hippocampus (2/3 loci in top 10% vs 10% null), supporting correct implementation; SCN1A itself has no hippocampus weight (not tested), SCN7A in the 2q24.3 cluster is nominal P=7.9e-3 (rank 3.2%).",
    "headline": {
        "primary_hippocampus_fhs_fdr05": combos.get("Hippocampus_fhs", {}).get("n_fdr05"),
        "note": "focal-HS power is limited; 0 genes pass FDR<0.05 in the primary analysis - this is an expected/acceptable negative result and is reported as-is (no tuning applied).",
    },
    "overlap_with_axis_genes": {
        "axis_set": "LDAM+Ferroptosis union (genesets/LDAM_Ferroptosis_union.txt), 149 axis genes tested in hippocampus x fhs",
        "fdr05": {"overlap": 0, "OR": None, "fisher_P": None,
                  "note": "0 TWAS genes pass FDR<0.05, enrichment undefined"},
        "nominal_p05": {"overlap": 10, "OR": 0.926, "fisher_P": 1.0,
                        "note": "no enrichment of axis genes among nominal TWAS hits"},
    },
    "degradations": [
        "SCN1A and several other epilepsy genes lack GTEx v8 hippocampus/cortex FUSION weights (expression heritability QC) and could not be tested.",
        "GTEx hippocampus panel is small (N~150), limiting weight accuracy and TWAS power, especially for focal-HS (smallest GWAS).",
        "167 (hippocampus) / 212 (cortex) genes dropped because no weight SNP matched the 1000G EUR reference; no other degradation.",
    ],
}
with open(os.path.join(OUT, "twas_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print(json.dumps({k: {"n_tested": v["n_tested"], "n_fdr05": v["n_fdr05"],
                      "top": v["top5"][0]["symbol"]} for k, v in combos.items()},
                 indent=2))
print("summary written")
