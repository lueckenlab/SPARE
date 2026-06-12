"""Compute a joint UMAP per CD45 compartment from the integrated embedding
(X_scVI_sample) of the *processed* h5ad and colour it by the important
metadata (cell type, response, timepoint, patient). Run after preprocess.

Usage: python make_umaps.py scct_db_cd45pos
"""
import os, sys, numpy as np, pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark"
ds = sys.argv[1] if len(sys.argv) > 1 else "scct_db_cd45pos"
d = os.path.join(REPO, "data", ds)
adata = sc.read_h5ad(os.path.join(d, f"{ds}_processed.h5ad"))

# attach sample-level metadata (response/timepoint/patient) by sample_id
meta = pd.read_csv(os.path.join(d, f"{ds}_metadata.csv"))
scol = "sample_id" if "sample_id" in meta.columns else meta.columns[0]
meta = meta.set_index(scol)
for c in ["response", "timepoint", "patient"]:
    if c in meta.columns:
        adata.obs[c] = adata.obs["sample_id"].map(meta[c]).astype("category")

# subsample for a fast, visually-identical UMAP (full 319k cells is slow on 2 CPUs)
MAXN = 80000
if adata.n_obs > MAXN:
    import numpy as _np
    rng = _np.random.RandomState(42)
    idx = rng.choice(adata.n_obs, MAXN, replace=False)
    adata = adata[idx].copy()
    print(f"subsampled to {adata.n_obs} cells for UMAP")

emb = "X_scVI_sample" if "X_scVI_sample" in adata.obsm else "X_pca"
print(ds, adata.shape, "| embedding:", emb)
sc.pp.neighbors(adata, use_rep=emb, n_neighbors=15)
sc.tl.umap(adata)

color_by = [c for c in ["cell_type", "response", "timepoint", "patient"] if c in adata.obs]
n = len(color_by)
fig, axes = plt.subplots(1, n, figsize=(6 * n, 5.2))
if n == 1:
    axes = [axes]
for ax, c in zip(axes, color_by):
    legend = "right margin" if adata.obs[c].nunique() <= 22 else None
    sc.pl.umap(adata, color=c, ax=ax, show=False, frameon=False,
               legend_loc=legend, size=4, title=c)
comp = "CD45+ (immune)" if ds.endswith("pos") else "CD45- (tumor/stroma)"
fig.suptitle(f"{ds} — {comp} — UMAP of X_scVI_sample", y=1.02, fontsize=13)
fig.tight_layout()
out = os.path.join(d, f"{ds}_umap.png")
fig.savefig(out, dpi=130, bbox_inches="tight")
print("SAVED", out)
