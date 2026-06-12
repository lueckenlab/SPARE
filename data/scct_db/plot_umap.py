"""Plot the precomputed scCT-DB cell-level UMAP for the largest dataset
(scCT_dataset_154-I, TNBC + pembrolizumab, GSE246613, ~278k cells), coloured by
the important metadata: cell type, sampling time (baseline vs 3 weeks) and
patient. UMAP coordinates ship in umap_info.txt (computed per dataset).
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
DS = "scCT_dataset_154-I"
U = pd.read_csv(os.path.join(HERE, "raw", DS, "umap_info.txt"), sep="\t")
U = U.rename(columns={"UMAP_1": "x", "UMAP_2": "y", "cell_type": "celltype",
                      "patient_id": "patient", "sampling_time": "timepoint"})
# downsample for legible scatter
rng = np.random.RandomState(42)
if len(U) > 120000:
    U = U.iloc[rng.choice(len(U), 120000, replace=False)]
print(DS, "cells plotted:", len(U), "| cell types:", U.celltype.nunique(),
      "| patients:", U.patient.nunique(), "| timepoints:", sorted(U.timepoint.unique()))

def scatter(ax, color_col, title, legend=True, max_leg=18):
    cats = U[color_col].astype(str)
    order = cats.value_counts().index.tolist()
    cmap = plt.get_cmap("tab20")
    colors = {c: cmap(i % 20) for i, c in enumerate(order)}
    for c in order:
        m = cats == c
        ax.scatter(U.x[m], U.y[m], s=1.2, linewidths=0, alpha=.6,
                   color=colors[c], rasterized=True)
    ax.set_title(title, fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("UMAP1", fontsize=8); ax.set_ylabel("UMAP2", fontsize=8)
    if legend:
        handles = [plt.Line2D([0], [0], marker="o", linestyle="", markersize=5,
                              color=colors[c]) for c in order[:max_leg]]
        labs = [c if len(c) < 30 else c[:28] + "…" for c in order[:max_leg]]
        ax.legend(handles, labs, frameon=False, fontsize=6, loc="center left",
                  bbox_to_anchor=(1.0, 0.5), handletextpad=0.2)

fig, axes = plt.subplots(1, 3, figsize=(20, 6))
scatter(axes[0], "celltype", f"{DS}: cell type", max_leg=20)
scatter(axes[1], "timepoint", "sampling time (baseline vs 3 weeks)")
scatter(axes[2], "patient", "patient", legend=False)
fig.suptitle("scCT-DB GSE246613 (TNBC + pembrolizumab) — precomputed cell-level UMAP",
             fontsize=13, y=1.02)
fig.tight_layout()
out = os.path.join(HERE, "gse246613_umap.png")
fig.savefig(out, dpi=130, bbox_inches="tight")
print("SAVED", out)
