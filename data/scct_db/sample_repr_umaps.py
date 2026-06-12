"""Plot the *sample-level* UMAP of each representation method (precomputed in
the aggregated h5ad as `<method>_umap`) coloured by sample-level metadata.
One figure per compartment: rows = metadata vars, cols = methods.

Usage: python sample_repr_umaps.py scct_db_cd45pos
"""
import os, sys, numpy as np, pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark"
ds = sys.argv[1] if len(sys.argv) > 1 else "scct_db_cd45pos"
a = ad.read_h5ad(os.path.join(REPO, "output_scct_db", f"{ds}.sample_representation.output.h5ad"))

# one representative variant per must-have method (fall back if absent)
WANT = [
    ("mrvi_epochs_100_umap", "mrvi"),
    ("sampleclr_scanvi_sample_id_t05_umap", "sampleclr"),
    ("sampleclr_supervised_scvi_sample_id_ba_umap", "sampleclr_sup"),
    ("mofa_pca_umap", "mofa"),
    ("gloscope_pca_umap", "gloscope"),
    ("cell_group_composition_clr_euclidean_umap", "composition(CLR)"),
    ("pascient_raw_counts_umap", "pascient"),
    ("pseudobulk_pca_euclidean_umap", "pseudobulk"),
]
methods = [(k, lab) for k, lab in WANT if k in a.obsm]
META = ["response", "timepoint", "sort", "dataset"]
META = [m for m in META if m in a.obs and a.obs[m].nunique() > 1]

palettes = {}
def colors_for(col):
    cats = list(pd.Series(a.obs[col].astype(str)).unique())
    cmap = plt.get_cmap("tab10")
    return {c: cmap(i % 10) for i, c in enumerate(sorted(cats))}

nr, nc = len(META), len(methods)
fig, axes = plt.subplots(nr, nc, figsize=(2.7 * nc, 2.7 * nr), squeeze=False)
for r, meta in enumerate(META):
    cmap = colors_for(meta)
    for c, (key, lab) in enumerate(methods):
        ax = axes[r][c]
        xy = a.obsm[key]
        vals = a.obs[meta].astype(str).values
        for cat in sorted(set(vals)):
            m = vals == cat
            ax.scatter(xy[m, 0], xy[m, 1], s=22, alpha=.85, linewidths=0,
                       color=cmap[cat], label=cat)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        if r == 0:
            ax.set_title(lab, fontsize=10)
        if c == 0:
            ax.set_ylabel(meta, fontsize=11, fontweight="bold")
    # one legend per row (rightmost)
    axes[r][-1].legend(frameon=False, fontsize=7, loc="center left",
                       bbox_to_anchor=(1.0, 0.5), markerscale=1.3)
comp = "CD45+ (immune)" if ds.endswith("pos") else "CD45- (tumor/stroma)"
fig.suptitle(f"{ds} — {comp} — sample-level representation UMAPs", y=1.005, fontsize=13)
fig.tight_layout()
out = os.path.join(REPO, "data", ds, f"{ds}_sample_repr_umaps.png")
fig.savefig(out, dpi=135, bbox_inches="tight")
print("SAVED", out, "| methods:", [m[1] for m in methods], "| meta:", META)
