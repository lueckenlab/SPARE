import scanpy as sc
import numpy as np
from patient_representation.pp import is_count_data

## VIASH START
par = {
    "input": "data/ticatlas.h5ad",
    "output": "data/ticatlas_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

print("Reading data")
adata = sc.read_h5ad(par["input"])

print("Excluding samples with non count data")

# Some datasets have TPMs instead of raw counts: https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/Tumor-Immune-Cell-Atlas/issues/10
NON_COUNT_DATASETS = ["breast", "lung1", "melanoma1"]  # "liver2" mentioned in the issue has raw count samples as well
NON_COUNT_SAMPLES = ["LIV2_6", "LIV2_7", "LIV2_8", "LIV2_9", "LIV2_10", "LIV2_11"]

NON_COUNT_DATA = adata.obs["source"].isin(NON_COUNT_DATASETS) | adata.obs["patient"].isin(NON_COUNT_SAMPLES)

print("Excluding ", adata[NON_COUNT_DATA, :].obs["patient"].nunique(), "samples and ", adata[NON_COUNT_DATA, :].n_obs, "cells with non count data")

adata = adata[~NON_COUNT_DATA, :].copy()

try:
    print("X contains count data:", is_count_data(adata.X))
except Exception as e:
    print("Error checking if adata.X contains count data:", e)

try:
    print("raw.X contains count data:", is_count_data(adata.raw.X))
except Exception as e:
    print("Error checking if raw.X contains count data:", e)

print("Copying raw counts to layers")
adata.X = adata.raw.X.copy()

try:
    print("X contains count data:", is_count_data(adata.X))
except Exception as e:
    print("Error checking if adata.X contains count data:", e)

adata.obs.loc[
    adata.obs["gender"].isin(["unknown", "NA"]),
    "gender"
] = np.nan

del adata.raw

adata.layers["X_raw_counts"] = adata.X.copy()
print("adata.obsm['X_raw_counts'].shape", adata.layers["X_raw_counts"].shape)

print("Filtering poor cells and genes")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
print(adata)

print("Saving output")
adata.write(par["output"], compression=par["output_compression"])
