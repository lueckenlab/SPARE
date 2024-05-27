import random

import scanpy as sc
import scvi
import pandas as pd
from patient_representation.pp import is_count_data
import numpy as np

random.seed(42)
np.random.seed(42)

## VIASH START
par = {
    "input": "data/combat.h5ad",
    "output": "data/combat_processed.h5ad",
    "cell_type_key":"Annotation_major_subset",
    "batch_covariates":["scRNASeq_sample_ID","Pool_ID"],
    "output_compression": "gzip",
}
## VIASH END

#combat
# CELL_TYPE_KEY = "Annotation_major_subset"
# SAMPLE_KEY = "scRNASeq_sample_ID"
# BATCH_KEY = "Pool_ID"

#onek1k
# SAMPLE_KEY = "donor_id"
# CELL_TYPE_KEY = "cell_type"
# BATCH_KEY = "pool_number"

#hlca
# SAMPLE_KEY = "donor_id"
# CELL_TYPE_KEY = "cell_type"
# BATCH_KEY = "dataset"

ADATA_PATH = par["input"]
CELL_TYPE_KEY = par["cell_type_key"]
BATCH_COVARIATES = par["batch_covariates"]
print("Reading data")
adata = sc.read_h5ad(ADATA_PATH)


##########################
#all to bedeleted, was just for testing
# print("X contains count data:", is_count_data(adata.X))
# # Copy raw counts to obsm
# adata.obsm["raw"] = adata.X.copy()
# adata.layers["raw"] = adata.X.copy()

#error for combat and onek1k1, only for hlca
####only for hlca -> because of raw.X
# print("raw.X contains count data:", is_count_data(adata.raw.X))
# Copy raw counts to obsm
# adata.obsm["raw"] = adata.raw.X.copy()
# adata.layers["raw"] = adata.raw.X.copy()
######################


# Find highly-variable genes
print("Subsetting HVG")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, flavor="seurat_v3",
                            n_top_genes=3000, layer="raw")
adata = adata[:, adata.var.highly_variable].copy()

####not batch_key given to calc HVG for COMBAT
##if batch_key is needed to subset to HVG -> change in batch_covariates args needed
#below is for ONEK1K
#sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000, layer="raw", batch_key=BATCH_KEY)

#below is for HLCA
#sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000, batch_key=BATCH_KEY, layer="raw")

# before subsetting to HVG for HLCA and ONEK1K
# print("Normalizing data")
# sc.pp.normalize_total(adata, target_sum=1e4)
# print("Log-transforming data")
# sc.pp.log1p(adata)

print("adata.shape", adata.shape)
print("adata.layers['raw'].shape", adata.layers["raw"].shape)

# Run PCA
# not being done for onek1k
print("Running PCA")
sc.tl.pca(adata)

#scvi, scanvi not for hlca ?
for batch_key in BATCH_COVARIATES:
    print(f"Obtain scVI and scanVI for batch key: {batch_key}")
    # Run scVI
    scvi.model.SCVI.setup_anndata(adata, layer="raw", batch_key=batch_key)
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    embedding_name_scVI = f"X_scVI_{batch_key}"
    adata.obsm[embedding_name_scVI] = vae.get_latent_representation()
    print(f"Embedding stored: {embedding_name_scVI}")

    # Run scANVI
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=adata,
        labels_key=CELL_TYPE_KEY,
        unlabeled_category="nan"
    )
    lvae.train(max_epochs=20, n_samples_per_label=100)
    embedding_name_scANVI = f"X_scANVI_{batch_key}"
    adata.obsm[embedding_name_scANVI] = lvae.get_latent_representation()
    print(f"Embedding stored: {embedding_name_scANVI}")

adata.write(par["output"], compression=par["output_compression"])