import random

import ehrapy as ep
import pandas as pd
import numpy as np
from pathlib import Path
random.seed(42)
np.random.seed(42)

## VIASH START
par = {
    "input": ["combat_pca.csv", "combat_random.csv"],
    "output": "data/combat/combat_representations.h5ad",
    "metadata_path":"data/combat/combat_metadata.csv",
    "obs_columns":["Source", "Outcome", "Death28", "Institute", "Pool_ID"],
    "cell_type_key": "Annotation_major_subset",
}
## VIASH END

print("Reading metadata")
metadata = pd.read_csv(par["metadata_path"], index_col=0)

obs_only_columns = par["obs_columns"]
if par["cell_type_key"]:
    obs_only_columns += metadata.columns[metadata.columns.str.startswith(par["cell_type_key"])].tolist()

print("Converting to AnnData")
meta_adata = ep.ad.df_to_anndata(
    metadata.drop(columns=obs_only_columns),  # Setting `columns_obs_only` causes a weird error
)
meta_adata.obs = metadata[obs_only_columns]

print("Encoding")
meta_adata = ep.pp.encode(meta_adata, autodetect=True)

print("Calculating QC metrics")
obs_metric, var_metrics = ep.pp.qc_metrics(meta_adata)

print("Imputing missing values")
ep.pp.knn_impute(meta_adata, n_neighbors=5)

print("Calculating PCA")
ep.pp.pca(meta_adata)

print("Calculating neighbors")
ep.pp.neighbors(
    meta_adata,
    use_rep="X_pca",
    key_added="ehrapy_neighbors",
)

print("Calculating Leiden")
ep.tl.leiden(meta_adata, key_added="ehrapy_leiden", neighbors_key="ehrapy_neighbors")

print("Calculating metadata-based UMAP")
ep.tl.umap(meta_adata, neighbors_key="ehrapy_neighbors")
meta_adata.obsm["ehrapy_umap"] = meta_adata.obsm["X_umap"]

samples = list(meta_adata.obs_names)
print(f"Samples in metadata: {len(samples)}")

for representation in par["input"]:
    print(f"Processing representation {representation}")
    representation_df = pd.read_csv(representation, index_col=0)

    # Remove samples that are not in the representation
    samples = [sample for sample in samples if sample in representation_df.index]
    representation_df = representation_df.loc[samples][samples]
    print(f"Current number of samples: {len(samples)}")
    if len(samples) < meta_adata.n_obs:
        meta_adata = meta_adata[samples, :]

    # Remove the file extension
    representation_name = Path(representation).stem

    # It makes more sense to put distances to the obsp, but pp.neighbors expects it to be in obsm
    # So we put full distances matrix in obsm, and obsp will contain only the distances for nearest neighbors 
    meta_adata.obsm[f"{representation_name}_distances"] = representation_df.loc[samples][samples]
    
    print("Putting neighbors to obsp")
    ep.pp.neighbors(
        meta_adata,
        use_rep=f"{representation_name}_distances",
        key_added=f"{representation_name}_neighbors",
        metric="precomputed"
    )

    print("Calculating Leiden")
    ep.tl.leiden(meta_adata, key_added=f"{representation_name}_leiden", neighbors_key=f"{representation_name}_neighbors")

    print("Calculating UMAP")
    ep.tl.umap(meta_adata, neighbors_key=f"{representation_name}_neighbors")
    meta_adata.obsm[f"{representation_name}_umap"] = meta_adata.obsm["X_umap"]

print("Resulting AnnData")
print(meta_adata)
meta_adata.write_h5ad(par["output"])
