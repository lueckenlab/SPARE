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
    "accessible_metadata_columns":["Age", "Sex", "BMI", "Hospitalstay", "SARSCoV2PCR", "TimeSinceOnset"],
    "cell_type_key": "Annotation_major_subset",
}
## VIASH END

print("Reading metadata")
metadata = pd.read_csv(par["metadata_path"], index_col=0)

obs_only_columns = metadata.columns.drop(par["accessible_metadata_columns"])

print("Converting to AnnData")
meta_adata = ep.ad.df_to_anndata(
    metadata[par["accessible_metadata_columns"]],  # Setting `columns_obs_only` causes a weird error
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

meta_adata.uns["sample_representations"] = ["ehrapy"]

for representation in par["input"]:
    print(f"Processing representation {representation}")
    representation_df = pd.read_csv(representation, index_col=0)

    nan_frac = representation_df.isna().to_numpy().mean()
    if nan_frac > 0.1:
        print(f"Skipping {representation}: {nan_frac:.1%} of values are NaN")
        continue

    # Remove samples that are not in the representation
    samples = [sample for sample in samples if sample in representation_df.index]
    representation_df = representation_df.loc[samples][samples]
    print(f"Current number of samples: {len(samples)}")
    if len(samples) < meta_adata.n_obs:
        meta_adata = meta_adata[samples, :]

    # Remove the file extension
    representation_name = Path(representation).stem

    # gloscope and similar KDE-based methods can return tiny negative distances
    # from floating-point error; pp.neighbors(metric="precomputed") rejects them.
    representation_df = representation_df.clip(lower=0)

    # Some methods (e.g. pilot_gm_vae) emit a small number of NaN entries when
    # the per-pair distance computation fails. The earlier 10% threshold has
    # already filtered out mostly-broken matrices, so the remaining NaNs are
    # rare; fill them with the matrix-wide max so sklearn's
    # KNeighborsTransformer (which rejects any NaN) treats those pairs as
    # "as far apart as possible".
    nan_mask = representation_df.isna()
    if nan_mask.to_numpy().any():
        fill_value = float(representation_df.fillna(0).to_numpy().max())
        representation_df = representation_df.fillna(fill_value)
        print(f"  Filled {int(nan_mask.to_numpy().sum())} NaN distances with max={fill_value:.4f}")

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

    meta_adata.uns["sample_representations"].append(representation_name)

print("Resulting AnnData")
print(meta_adata)
meta_adata.write_h5ad(par["output"])
