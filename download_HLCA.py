import cellxgene_census
import scanpy as sc
import pandas as pd

FULL_HLCA_ID = "9f222629-9e39-47d0-b83f-e08d610c7479"

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

try:
    cellxgene_census.download_source_h5ad(
        FULL_HLCA_ID, to_path="../data/HLCA.h5ad"
    )
except Exception as e:
    print("Failed to download:", e)

adata = sc.read_h5ad("../data/HLCA.h5ad")

print("HLCA shape before filtering out diseases:", adata.shape)
adata = adata[adata.obs["disease"].isin(pneumonia_related_diseases)]
print("HLCA shape after filtering out diseases:", adata.shape)

print("Filtering cells with no annotation")
adata = adata[pd.notna(adata.obs["study"])]
print("HLCA shape after filtering out unannotated cells:", adata.shape)

adata.write("../data/HLCA_subset.h5ad")
