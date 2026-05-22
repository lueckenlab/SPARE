import random

import scanpy as sc
import pandas as pd
import patpy as pr
import numpy as np

random.seed(42)
np.random.seed(42)

## VIASH START
par = {
    "input": "data/combat/combat_processed.h5ad",
    "output": "data/combat_pilot_gm_vae_pca.csv",
    "sample_key": "scRNASeq_sample_ID",
    "cell_type_key": "Annotation_major_subset",
    "sample_state_col": None,
    "layer": "X_pca",
    "num_classes": 11,
    "gaussian_size": 64,
    "epochs": 50,
    "cuda": 0,
}
## VIASH END

print("Reading adata")
adata = sc.read(par["input"])

print(adata)

# PILOT-GM-VAE requires a sample "status" column. Internally pilotgm concatenates
# [PCA, sample_col, component_assignment, status] and renames the last three columns
# by position to sampleID/cell_type/status. If the status column shares a name with the
# sample column the rename collides and the 'sampleID' column is lost. The status only
# feeds optional ARI/silhouette evaluation (not the distances), so when no distinct
# status column is supplied we add a constant one with a dedicated name.
sample_state_col = par["sample_state_col"]
if not sample_state_col or sample_state_col == par["sample_key"]:
    sample_state_col = "_pilot_gm_vae_status"
    adata.obs[sample_state_col] = "all"

print("Setting up the representation method")
representation_method = pr.tl.PILOTGMVAE(
    sample_key=par["sample_key"],
    sample_state_col=sample_state_col,
    layer=par["layer"],
    num_classes=par["num_classes"],
    gaussian_size=par["gaussian_size"],
    epochs=par["epochs"],
    cuda=par["cuda"],
)

print("Preparing the anndata object (training GM-VAE)")
representation_method.prepare_anndata(adata)

print("Calculating the distance matrix")
distances = representation_method.calculate_distance_matrix(
    force=True
)

print("Saving results")
distances_df = pd.DataFrame(distances, index=representation_method.samples, columns=representation_method.samples)
distances_df.to_csv(par["output"], index=True)
