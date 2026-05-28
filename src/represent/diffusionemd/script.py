import random

import scanpy as sc
import pandas as pd
import patpy as pr
import numpy as np

# DiffusionEMD 0.5.0 imports `DistanceMetric` from sklearn.neighbors, which moved
# to sklearn.metrics in sklearn >= 1.3. Shim it so the (unmaintained) package imports.
import sklearn.neighbors
import sklearn.metrics
if not hasattr(sklearn.neighbors, "DistanceMetric"):
    sklearn.neighbors.DistanceMetric = sklearn.metrics.DistanceMetric

random.seed(42)
np.random.seed(42)

## VIASH START
par = {
    "input": "data/combat/combat_processed.h5ad",
    "output": "data/combat_diffusionemd_pca.csv",
    "sample_key": "scRNASeq_sample_ID",
    "cell_type_key": "Annotation_major_subset",
    "layer": "X_pca",
    "n_neighbors": 15,
    "n_scales": 6,
    "subset_n_cells_per_sample": 2000,
}
## VIASH END

print("Reading adata")
adata = sc.read(par["input"])

print(adata)

# Subsample to at most N cells per donor. The kNN graph scales with total
# n_cells, so this caps the graph size at n_donors * N. Donors with fewer
# cells than the cap are kept whole.
n_per_sample = par["subset_n_cells_per_sample"]
if n_per_sample and n_per_sample > 0:
    print(f"Subsampling to ≤{n_per_sample} cells per donor")
    rng = np.random.default_rng(42)
    sample_ids = adata.obs[par["sample_key"]]
    keep_idxs = []
    for level in sample_ids.unique():
        idxs = np.where((sample_ids == level).to_numpy())[0]
        if len(idxs) <= n_per_sample:
            keep_idxs.extend(idxs.tolist())
        else:
            keep_idxs.extend(rng.choice(idxs, size=n_per_sample, replace=False).tolist())
    adata = adata[np.sort(np.asarray(keep_idxs))].copy()
    print(f"After subsampling: {adata.n_obs:,d} cells across {adata.obs[par['sample_key']].nunique()} donors")

print("Setting up the representation method")
representation_method = pr.tl.DiffusionEarthMoverDistance(
    sample_key=par["sample_key"],
    cell_group_key=par["cell_type_key"],
    layer=par["layer"],
    n_neighbors=par["n_neighbors"],
    n_scales=par["n_scales"],
)

print("Preparing the anndata object")
representation_method.prepare_anndata(adata)

print("Calculating the distance matrix")
distances = representation_method.calculate_distance_matrix(
    force=True
)

print("Saving results")
distances_df = pd.DataFrame(distances, index=representation_method.samples, columns=representation_method.samples)
distances_df.to_csv(par["output"], index=True)
