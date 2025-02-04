import scanpy as sc
import pandas as pd
from patient_representation.pp import is_count_data

## VIASH START
par = {
    "input": "data/hlca.h5ad",
    "output": "data/hlca_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

# The whole list of diseases in HLCA. The ones that are not pneumonia-related are commented out.
pneumonia_related_diseases = [
    "disease",
    "normal",
    "COVID-19",
    "pulmonary fibrosis"
    "interstitial lung disease",
    # chronic obstructive pulmonary disease
    # lung adenocarcinoma
    "pneumonia",
    # chronic rhinitis
    # lung large cell carcinoma
    # squamous cell lung carcinoma
    # cystic fibrosis
    # lymphangioleiomyomatosis
    # pleomorphic carcinoma
    # hypersensitivity pneumonitis
    "non-specific interstitial pneumonia",
    # pulmonary sarcoidosis
]

adata = sc.read_h5ad(par["input"])
adata.X = adata.raw.X.copy()

print("HLCA shape before filtering out diseases:", adata.shape)
adata = adata[adata.obs["disease"].isin(pneumonia_related_diseases)]
print("HLCA shape after filtering out diseases:", adata.shape)

print("Filtering cells with no annotation")
adata = adata[pd.notna(adata.obs["study"])]
print("HLCA shape after filtering out unannotated cells:", adata.shape)

print("X contains count data:", is_count_data(adata.X))
print("raw.X contains count data:", is_count_data(adata.raw.X))

# Copy counts to layers
adata.layers["X_raw_counts"] = adata.X.copy()
del adata.raw
print(adata)

adata.write(par["output"], compression=par["output_compression"])
