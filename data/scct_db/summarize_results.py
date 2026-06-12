"""Summarize response-prediction results across both CD45 compartments.
Produces a merged per-method `relevant` (response) score table + a comparison
figure, highlighting the must-have methods. Run after evaluate.
"""
import os, re, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
EVAL = os.path.join(REPO, "output_scct_db_evaluate")
MUST = ["mrvi", "cell_group_composition", "sampleclr", "mofa", "gloscope",
        "pascient", "sampleclr_supervised"]

def method_of(rep):
    for m in sorted(MUST, key=len, reverse=True):
        if rep.startswith(m):
            return m
    return rep.split("_")[0]

frames = {}
for ds in ["scct_db_cd45pos", "scct_db_cd45neg"]:
    d = pd.read_csv(os.path.join(EVAL, ds, "averaged_scores.csv"))
    d = d[["representation", "relevant", "technical", "total"]].copy()
    d["compartment"] = "CD45+" if ds.endswith("pos") else "CD45-"
    frames[ds] = d

alld = pd.concat(frames.values(), ignore_index=True)
alld["method"] = alld["representation"].map(method_of)

# Best variant per method per compartment, on the response (relevant) axis
best = (alld.sort_values("relevant", ascending=False)
        .groupby(["compartment", "method"], as_index=False).first())
piv = best.pivot(index="method", columns="compartment",
                 values="relevant").sort_values("CD45+", ascending=False)
piv.to_csv(os.path.join(HERE, "response_scores_by_method.csv"))
print("=== Best response (relevant) score per method ===")
print(piv.round(3).to_string())

# Full ranking table (top variants) per compartment
for ds, d in frames.items():
    top = d.sort_values("relevant", ascending=False).head(15)
    top.to_csv(os.path.join(HERE, f"{ds}_response_ranking.csv"), index=False)

# ---- figure: response score per method, both compartments ----
order = piv.index.tolist()
fig, ax = plt.subplots(figsize=(9, 6))
y = np.arange(len(order)); h = 0.38
for i, (comp, color) in enumerate([("CD45+", "#4C72B0"), ("CD45-", "#C44E52")]):
    vals = [piv.loc[m, comp] if comp in piv.columns and not pd.isna(piv.loc[m, comp]) else 0
            for m in order]
    ax.barh(y + (0.5 - i) * h, vals, height=h, label=comp, color=color)
# random-vector floor reference
floor = alld[alld.representation.str.startswith("random_vector")]["relevant"].max()
ax.axvline(floor, ls="--", color="grey", lw=1)
ax.text(floor, len(order) - 0.3, f" random floor {floor:.2f}", color="grey", fontsize=8, va="top")
ax.set_yticks(y); ax.set_yticklabels(order)
ax.invert_yaxis()
ax.set_xlim(0, 1); ax.set_xlabel("response score (relevant, best variant)")
ax.set_title("scCT-DB GSE246613 — response prediction by method")
ax.legend(frameon=False, title="compartment")
sns.despine(ax=ax); ax.grid(False)
fig.tight_layout()
out = os.path.join(HERE, "response_prediction_summary.png")
fig.savefig(out, dpi=140)
print("SAVED", out)
