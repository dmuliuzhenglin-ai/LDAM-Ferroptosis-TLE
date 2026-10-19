import pandas as pd, os, ast, pickle
from ctxcore.genesig import Regulon

BASE = os.getcwd()
OUT = os.path.join(BASE, "results", "56_pySCENIC")

def main():
    df = pd.read_csv(os.path.join(OUT, "regulons_df.tsv"), sep="\t",
                     header=[0, 1], index_col=[0, 1])
    # flatten columns to field names (second header level)
    df.columns = [c[1] for c in df.columns]
    regs = []
    for tf_name, grp in df.groupby(level=0):
        targets = {}
        activating = True
        nes_max = 0.0
        for _, r in grp.iterrows():
            ctx_str = str(r["Context"])
            if "repressing" in ctx_str:
                activating = False
            try:
                nes_max = max(nes_max, float(r["NES"]))
            except Exception:
                pass
            tg = r["TargetGenes"]
            if isinstance(tg, str):
                tg = ast.literal_eval(tg)
            for gene, w in tg:
                if gene not in targets or w > targets[gene]:
                    targets[gene] = w
        if len(targets) >= 5:
            regs.append(Regulon(
                name=f"{tf_name}({'+' if activating else '-'})",
                gene2weight=targets,
                gene2occurrence=[],
                transcription_factor=str(tf_name),
                context=frozenset({"activating" if activating else "repressing"}),
                nes=nes_max))
    with open(os.path.join(OUT, "regulons.pkl"), "wb") as f:
        pickle.dump(regs, f)
    print("regulons:", len(regs))
    print("sample:", regs[0].name, len(regs[0].genes), list(regs[0].genes)[:5])

if __name__ == "__main__":
    main()
