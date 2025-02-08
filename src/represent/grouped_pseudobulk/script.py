import random

import scanpy as sc
import pandas as pd
import patient_representation as pr
import numpy as np

random.seed(42)
np.random.seed(42)

## VIASH START
par = {
    "input": "data/combat_processed.h5ad",
    "output": "data/combat_pseudobulk_pca.h5ad",
    "sample_key":"scRNASeq_sample_ID",
    "cell_type_key": "Annotation_major_subset",
    "layer": "X_pca",
    "aggregation": "mean",
    "distance_metric": "euclidean",
}
## VIASH END

print("Reading adata")
adata = sc.read(par["input"])

print(adata)

print("Setting up th representation method")
representation_method = pr.tl.GroupedPseudobulk(
    sample_key=par["sample_key"],
    cell_group_key=par["cell_type_key"],
    layer=par["layer"]
)

print("Preparing the anndata object")
representation_method.prepare_anndata(adata)

print("Calculating the distance matrix")
distances = representation_method.calculate_distance_matrix(
    force=True,
    aggregate=par["aggregation"],
    dist=par["distance_metric"]
)

print("Saving results")
distances_df = pd.DataFrame(distances, index=representation_method.samples, columns=representation_method.samples)
distances_df.to_csv(par["output"], index=True)
