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
}
## VIASH END

print("Reading adata")
adata = sc.read(par["input"])

print(adata)

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
