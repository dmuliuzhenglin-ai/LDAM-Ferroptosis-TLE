"""Phase 2b: MR + coloc for LPCAT2, CD9, S100A10 (MetaBrain hippocampus) vs ILAE GWAS."""
import sys, os, gzip, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from mr_coloc_utils import ivw, mr_egger, wald_ratio, coloc_abf, clump
from p2_lgals3 import load_gwas, harmonize, run_mr, pick_instruments  # reuse

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
RES = os.path.join(BASE, "results", "42_MR_coloc")

genes = {"LPCAT2": 16, "CD9": 12, "S100A10": 1}
all_res = []
for gene in genes:
    f = os.path.join(RES, f"exposure_{gene}_MetaBrain_hippocampus_raw.tsv")
    df = pd.read_csv(f, sep="\t")
    rows = []
    for _, r in df.iterrows():
        snp = r["SNP"].split(":")
        if len(snp) >= 3 and snp[2].startswith("rs"):
            rows.append(dict(SNP=snp[2], chr=int(r["SNPChr"]), pos=int(r["SNPPos"]),
                             EA=r["SNPEffectAllele"], eaf=float(r["SNPEffectAlleleFreq"]),
                             beta=float(r["MetaBeta"]), se=float(r["MetaSE"]), p=float(r["MetaP"])))
    ex = pd.DataFrame(rows).drop_duplicates(subset=["SNP"])
    ex.to_csv(os.path.join(RES, f"exposure_{gene}_MetaBrain_hippocampus.csv"), index=False)
    gw = load_gwas(ex["SNP"].unique())
    m = harmonize(ex, gw, "EA", "eaf")
    m.to_csv(os.path.join(RES, f"harmonized_{gene}_hippocampus_vs_ILAEall.csv"), index=False)
    inst, mode = pick_instruments(m, "p", "pos", "chr")
    mr = run_mr(inst, "beta", "se") if len(inst) else {}
    co = coloc_abf(m["beta"].values, m["se"].values, m["beta_gwas_h"].values, m["se_gwas"].values) if len(m) else {}
    print(gene, "| harm:", len(m), "| instr:", len(inst), mode)
    print("  MR:", json.dumps(mr, default=float)[:400])
    print("  coloc:", json.dumps(co, default=float)[:300])
    all_res.append(dict(gene=gene, tissue="MetaBrain hippocampus EUR (N=168)",
                        n_harmonized=len(m), n_instruments=len(inst), instrument_mode=mode,
                        **{f"MR:{k}": v for k, v in mr.items()}, coloc=co))

with open(os.path.join(RES, "hub_genes_results.json"), "w") as f:
    json.dump(all_res, f, indent=1, default=float)
