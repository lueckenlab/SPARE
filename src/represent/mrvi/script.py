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
    "output": "data/combat_mrvi.csv",
    "cell_type_key":"Annotation_major_subset",
    "sample_key":"scRNASeq_sample_ID",
    "max_epochs": 100,
    "use_gpu": True,
}
## VIASH END

print("Reading adata")
adata = sc.read(par["input"])
print(adata)

print("Setting up the representation method")
representation_method = pr.tl.MrVI(
    sample_key=par["sample_key"],
    cell_group_key=par["cell_type_key"],
    layer="X_raw_counts",
    max_epochs=par["max_epochs"],
)

print("Preparing the anndata object")
representation_method.prepare_anndata(adata)

print("Calculating the distance matrix")
distances = representation_method.calculate_distance_matrix(force=True)

print("Saving results")
distances_df = pd.DataFrame(distances, index=representation_method.samples, columns=representation_method.samples)
distances_df.to_csv(par["output"], index=True)
