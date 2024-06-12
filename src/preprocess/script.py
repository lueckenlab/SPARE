import random

import scanpy as sc
import scvi
import pandas as pd
from patient_representation.pp import is_count_data
import numpy as np
import os
import matplotlib.pyplot as plt

random.seed(42)
np.random.seed(42)

## VIASH START
par = {
    "input": "data/combat.h5ad",
    "output": "data/combat_processed.h5ad",
    "cell_type_key":"Annotation_major_subset",
    "batch_covariates":["scRNASeq_sample_ID","Pool_ID"],
    "batch_effect_covariate":"Pool_ID",
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

#stephenson
# SAMPLE_KEY = "sample_id"
# CELL_TYPE_KEY = "cell_type"
# BATCH_KEY = "Site"

ADATA_PATH = par["input"]
CELL_TYPE_KEY = par["cell_type_key"]
BATCH_COVARIATES = par["batch_covariates"]
BATCH_EFFECT = par["batch_effect_covariate"]

figures_directory = os.path.join(os.path.dirname(par["output"]), "figures")

if not os.path.exists(figures_directory):
    os.makedirs(figures_directory)

print("Reading data")
adata = sc.read_h5ad(ADATA_PATH)

print("Normalizing data")
sc.pp.normalize_total(adata, target_sum=1e4)
print("Log-transforming data")
sc.pp.log1p(adata)

# Find highly-variable genes
print("Subsetting HVG")
sc.pp.highly_variable_genes(adata, flavor="seurat_v3",span=0.5, n_top_genes=3000, layer="raw", batch_key=BATCH_EFFECT)
adata = adata[:, adata.var.highly_variable].copy()


print("adata.shape", adata.shape)
print("adata.layers['raw'].shape", adata.layers["raw"].shape)

print("Running PCA")
sc.tl.pca(adata)

def plot_loss(history, loss_keys, title, filenames, counter):
    plt.figure(figsize=(10, 5))
    for loss_key, filename in zip(loss_keys, filenames):
        if loss_key in history:
            plt.plot(history[loss_key], label=f'{loss_key} (Train)')
            plt.xlabel('Epoch')
            plt.ylabel('Loss')
            plt.title(title)
            plt.legend()
            # Add counter to filename
            filename_with_counter = f"{counter}_{filename}"
            full_path = os.path.join(figures_directory, filename_with_counter)
            plt.savefig(full_path)
            plt.close()
            counter += 1
    return counter

counter = 1
for batch_key in BATCH_COVARIATES:
    print(f"Obtain scVI and scanVI for batch key: {batch_key}")
    # Run scVI
    scvi.model.SCVI.setup_anndata(adata, layer="raw", batch_key=batch_key)
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    embedding_name_scVI = f"X_scVI_{batch_key}"
    adata.obsm[embedding_name_scVI] = vae.get_latent_representation()
    print(f"Embedding stored: {embedding_name_scVI}")

    counter = plot_loss(vae.history, 
        ['reconstruction_loss_train', 'train_loss_epoch', 'elbo_train'], 
        f'Training Loss Metrics for scVI ({batch_key})', 
        [f'reconstruction_loss_scVI_{batch_key}.png', f'train_loss_epoch_scVI_{batch_key}.png', f'elbo_train_scVI_{batch_key}.png'], counter)

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

    counter = plot_loss(vae.history, 
    ['reconstruction_loss_train', 'train_loss_epoch', 'elbo_train'], 
    f'Training Loss Metrics for scANVI ({batch_key})', 
    [f'reconstruction_loss_scANVI_{batch_key}.png', f'train_loss_epoch_scANVI_{batch_key}.png', f'elbo_train_scANVI_{batch_key}.png'], counter)


adata.write(par["output"], compression=par["output_compression"])
