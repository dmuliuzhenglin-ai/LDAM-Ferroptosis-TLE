"""B4 figure: forest plot of hub-gene pQTL/eQTL MR estimates vs epilepsy."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"C:\Users\19367\Desktop\小胶质细胞脂滴累积（LDAM）联合铁死亡在海马起源颞叶癫痫中的生信研究"
FIGD = os.path.join(BASE, "figures", "replication_v24")
d = pd.read_csv(os.path.join(BASE, "results", "54_hub_pQTL_MR", "hub_pqtl_mr_results.csv"))
d = d[d["status"] == "ok"].copy()
d["lo"] = d["beta"] - 1.96 * d["se"]; d["hi"] = d["beta"] + 1.96 * d["se"]
d["lab"] = d["gene"] + " (" + d["layer"].str.replace("pQTL_deCODE", "pQTL").str.replace("eQTL_eQTLGen_blood", "eQTL") + ", n=" + d["n_snp"].astype(int).astype(str) + ")"

fig, ax = plt.subplots(figsize=(7, 4.2))
order = d.sort_values(["outcome", "gene"], ascending=[False, True]).reset_index(drop=True)
y = np.arange(len(order))[::-1]
colors = ["#2c7fb0" if o == "all_epilepsy_ILAE2018" else "#cb181d" for o in order["outcome"]]
ax.errorbar(order["beta"], y, xerr=[order["beta"] - order["lo"], order["hi"] - order["beta"]],
            fmt="none", ecolor="#bbbbbb", zorder=1)
ax.scatter(order["beta"], y, c=colors, s=40, zorder=2)
ax.axvline(0, color="k", ls="--", lw=0.7)
ax.set_yticks(y)
ax.set_yticklabels(order["lab"] + "  |  " + order["outcome"].str.replace("_ILAE2018", "").str.replace("_", " "), fontsize=7)
ax.set_xlabel("MR effect on epilepsy risk (beta, log-OR scale)")
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(color="#2c7fb0", label="All epilepsy"),
                   mpatches.Patch(color="#cb181d", label="Focal epilepsy")], fontsize=7, frameon=False)
ax.set_title("Hub-gene pQTL/eQTL Mendelian randomization", fontsize=9)
plt.tight_layout()
fig.savefig(os.path.join(FIGD, "Figure_hubMR_v24.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIGD, "Figure_hubMR_v24.pdf"), bbox_inches="tight")
print("saved")
