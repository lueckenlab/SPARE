import random

import scanpy as sc
import pandas as pd
import patpy as pr
import numpy as np

# PhEMD is defined in patpy but not re-exported under patpy.tl, so import it directly
from patpy.tl.sample_representation import PhEMD

random.seed(42)
np.random.seed(42)

## VIASH START
par = {
    "input": "data/combat/combat_processed.h5ad",
    "output": "data/combat_phemd_pca.csv",
    "sample_key": "scRNASeq_sample_ID",
    "cell_type_key": "Annotation_major_subset",
    "layer": "X_pca",
    "n_clusters": 8,
    "subset_fraction": None,
    "subset_n_obs": None,
    "subset_min_obs_per_sample": 500,
    "n_jobs": -1,
}
## VIASH END

print("Reading adata")
adata = sc.read(par["input"])

print(adata)

# PHATE is slow, so optionally subsample cells per cell type first. We do this here
# rather than via PhEMD.prepare_anndata(subset_*=...) / patpy.pp.subsample because that
# path is doubly broken in the pinned patpy: a stale kwarg name (min_obs_per_category)
# and no clamping of n_obs to the available cells in a category (it samples without
# replacement and crashes when a cell type has fewer than n_obs cells). Cell types with
# <= subset_min_obs_per_sample cells are kept whole; otherwise we draw min(target, available).
if par["subset_fraction"] is not None or par["subset_n_obs"] is not None:
    print("Subsampling cells per cell type before PhEMD")
    rng = np.random.default_rng(42)
    cell_types = adata.obs[par["cell_type_key"]]
    floor = par["subset_min_obs_per_sample"]
    keep_idxs = []
    for level in cell_types.unique():
        idxs = np.where((cell_types == level).to_numpy())[0]
        n_available = len(idxs)
        if n_available <= floor:
            keep_idxs.extend(idxs.tolist())
            continue
        if par["subset_fraction"] is not None:
            target = int(par["subset_fraction"] * n_available)
        else:
            target = int(par["subset_n_obs"])
        target = max(1, min(target, n_available))
        keep_idxs.extend(rng.choice(idxs, size=target, replace=False).tolist())
    adata = adata[np.sort(np.asarray(keep_idxs))].copy()
    print(f"Subsampled to {adata.n_obs} cells")

print("Setting up the representation method")
representation_method = PhEMD(
    sample_key=par["sample_key"],
    cell_group_key=par["cell_type_key"],
    layer=par["layer"],
    n_clusters=par["n_clusters"],
)

print("Preparing the anndata object")
representation_method.prepare_anndata(adata)

print("Calculating the distance matrix")
distances = representation_method.calculate_distance_matrix(
    force=True,
    n_jobs=par["n_jobs"],
)

print("Saving results")
distances_df = pd.DataFrame(distances, index=representation_method.samples, columns=representation_method.samples)
distances_df.to_csv(par["output"], index=True)
