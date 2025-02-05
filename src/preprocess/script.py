import random

import scanpy as sc
import scvi
import pandas as pd
import patient_representation as pr
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
    "sample_key":"scRNASeq_sample_ID",
    "batch_covariates":["scRNASeq_sample_ID","Pool_ID","Institute"],
    "samples_metadata_cols":["Source", "Outcome", "Death28", "Institute", "Pool_ID"],
    "batch_key":"Pool_ID",
    "output_compression": "gzip",
    "sample_size_threshold": 100,
    "output_metadata": "data/combat_metadata_200.csv",
}
## VIASH END

CELL_TYPE_KEY = par["cell_type_key"]
BATCH_KEY = par["batch_key"]  # Used for HVG selection and integration
SAMPLE_KEY = par["sample_key"]

figures_directory = os.path.join(os.path.dirname(par["output"]), "figures")

if not os.path.exists(figures_directory):
    os.makedirs(figures_directory)

print("Reading data")
adata = sc.read_h5ad(par["input"])

print("Filtering small samples")
adata = pr.pp.filter_small_samples(adata, sample_key=SAMPLE_KEY, sample_size_threshold=par["sample_size_threshold"])

print("Normalizing data")
sc.pp.normalize_total(adata, target_sum=1e4)
print("Log-transforming data")
sc.pp.log1p(adata)

# Find highly-variable genes
print("Subsetting HVG")
sc.pp.highly_variable_genes(adata, flavor="seurat_v3",span=0.5, n_top_genes=3000, layer="X_raw_counts", batch_key=BATCH_KEY)
adata = adata[:, adata.var.highly_variable].copy()


print("adata.shape", adata.shape)
print("adata.layers['X_raw_counts'].shape", adata.layers["X_raw_counts"].shape)

print("Running PCA")
sc.tl.pca(adata)

def plot_loss(history, loss_keys, title, filenames, counter):
    plt.figure(figsize=(10, 5))
    for loss_key, filename in zip(loss_keys, filenames):
        if loss_key in history:
            plt.plot(history[loss_key], label=f"{loss_key} (Train)")
            plt.xlabel("Epoch")
            plt.ylabel("Loss")
            plt.title(title)
            plt.legend()
            filename_with_counter = f"{counter}_{filename}"
            full_path = os.path.join(figures_directory, filename_with_counter)
            plt.savefig(full_path)
            plt.close()
            counter += 1
    return counter

counter = 1
for batch_key in par["batch_covariates"]:
    # print("Running Harmony on PCA embeddings")
    # sc.external.pp.harmony_integrate(adata, basis="X_pca", key=batch_key, adjusted_basis=f"X_harmony_{batch_key}")

    print(f"Obtain scVI and scanVI for batch key: {batch_key}")
    # Run scVI
    scvi.model.SCVI.setup_anndata(adata, layer="X_raw_counts", batch_key=batch_key)
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    embedding_name_scVI = f"X_scVI_{batch_key}"
    adata.obsm[embedding_name_scVI] = vae.get_latent_representation()
    print(f"Embedding stored: {embedding_name_scVI}")

    counter = plot_loss(
        history=vae.history, 
        loss_keys=["reconstruction_loss_train", "train_loss_epoch", "elbo_train"], 
        title=f"Training Loss Metrics for scVI ({batch_key})", 
        filenames=[f"reconstruction_loss_scVI_{batch_key}.png", f"train_loss_epoch_scVI_{batch_key}.png", f"elbo_train_scVI_{batch_key}.png"],
        counter=counter
    )

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

    counter = plot_loss(
        history=lvae.history, 
        loss_keys=["reconstruction_loss_train", "train_loss_epoch", "elbo_train"], 
        title=f"Training Loss Metrics for scANVI ({batch_key})", 
        filenames=[f"reconstruction_loss_scANVI_{batch_key}.png", f"train_loss_epoch_scANVI_{batch_key}.png", f"elbo_train_scANVI_{batch_key}.png"],
        counter=counter
    )

print("Running scPoli")
scpoli = pr.tl.SCPoli(sample_key=SAMPLE_KEY, cell_group_key=CELL_TYPE_KEY, layer="X_raw_counts")

print("Preparing adata")
scpoli.prepare_anndata(adata, optimize_adata=False)

print("Calculating distances")
scpoli_distances = scpoli.calculate_distance_matrix(force=True)

print("Saving scPoli distances")
adata.uns["scpoli_distances"] = scpoli_distances
adata.uns["scpolis_samples"] = scpoli.samples

print("Saving scPoli cell representation")
adata.obsm["X_scpoli"] = scpoli.model.get_latent(
    scpoli.adata,
    mean=True
)

# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True
)

cell_qc_metadata = pr.pp.calculate_cell_qc_metrics(adata, sample_key=SAMPLE_KEY, cell_qc_vars=["n_genes_by_counts", "total_counts", "pct_counts_mt"])
n_genes_metadata = pr.pp.calculate_n_cells_per_sample(adata, SAMPLE_KEY)
composition_metadata = pr.pp.calculate_compositional_metrics(adata, SAMPLE_KEY, [CELL_TYPE_KEY], normalize_to=100)

metadata = pr.pp.extract_metadata(adata, SAMPLE_KEY, par["samples_metadata_cols"])
metadata = pd.concat([
    metadata,
    cell_qc_metadata.loc[metadata.index],
    n_genes_metadata.loc[metadata.index],
    composition_metadata.loc[metadata.index]
], axis=1)

print("ADATA PREPROCESSED: ")
print(adata)

adata.write(par["output"], compression=par["output_compression"])
metadata.to_csv(par["output_metadata"])
