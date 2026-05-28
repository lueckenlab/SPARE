import random

import scanpy as sc
import pandas as pd
import patpy as pr
import numpy as np

random.seed(42)
np.random.seed(42)

## VIASH START
par = {
    "input": "data/combat_processed.h5ad",
    "output": "data/combat_gloscope_pca.csv",
    "sample_key":"scRNASeq_sample_ID",
    "cell_type_key": "Annotation_major_subset",
    "layer": "X_pca",
    "use_gpu": True,
    "n_components": 10,
}
## VIASH END

print("Reading adata")
adata = sc.read(par["input"])

print(adata)

print("Setting up the representation method")
representation_method = pr.tl.GloScope_py(
    sample_key=par["sample_key"],
    cell_group_key=par["cell_type_key"],
    layer=par["layer"],
    use_gpu=par["use_gpu"],
    n_components=par.get("n_components"),
    seed=42,
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
