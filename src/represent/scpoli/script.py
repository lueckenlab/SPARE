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
    "max_epochs": 50,
}
## VIASH END

print("Reading adata")
adata = sc.read(par["input"])

print("Running scPoli")
if "scpoli_distances" in adata.uns and "scpolis_samples" in adata.uns:
    print("Using existing scPoli distances")
    distances = adata.uns["scpoli_distances"]
    samples = adata.uns["scpolis_samples"]
else:
    print("Running scPoli")
    representation_method = pr.tl.SCPoli(
        sample_key=par["sample_key"],
        cell_group_key=par["cell_type_key"],
        layer="X_raw_counts",
        n_epochs=par["max_epochs"],
    )
    print("Preparing adata")
    representation_method.prepare_anndata(adata)

    distances = representation_method.calculate_distance_matrix(force=True)
    samples = representation_method.samples

distances_df = pd.DataFrame(distances, index=samples, columns=samples)
distances_df.to_csv(par["output"], index=True)
