import cellxgene_census

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

with cellxgene_census.open_soma() as census:
    adata = cellxgene_census.get_anndata(
        census, organism="Homo sapiens", obs_value_filter=f"dataset_id == '{FULL_HLCA_ID}' and disease in {str(pneumonia_related_diseases)}"
    )

adata.write("data/HLCA_subset.h5ad")