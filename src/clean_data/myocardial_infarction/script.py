"""Clean the Kuppe et al. 2022 human myocardial infarction snRNA-seq atlas.

The file ships from CELLxGENE (schema 7.0.0): normalized counts in `.X`,
raw integer counts in `.raw`, genes indexed by Ensembl IDs with symbols in
`var["feature_name"]`. This mirrors the lupus / sound_life cleaners:

1. Rebuild an AnnData with raw counts in `.X` (src/preprocess re-normalizes
   from there and creates the X_raw_counts / X_log_norm_counts layers).
2. Use gene symbols as var_names for consistency with the other datasets.
3. Keep a curated, sample-level-meaningful obs; derive a numeric `age` and
   an ordinal `disease_progression` (myogenic 0 -> ischemic 1 -> fibrotic 2)
   used as the trajectory variable in evaluate.
"""
import re

import anndata as ad
import numpy as np
import scanpy as sc

## VIASH START
par = {
    "input": "data/myocardial_infarction/myocardial_infarction.h5ad",
    "output": "data/myocardial_infarction/myocardial_infarction_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

# Columns worth carrying forward (the rest are per-cell QC scores or
# ontology-term-id duplicates of the human-readable columns).
KEEP_OBS = [
    "sample",
    "donor_id",
    "patient_region_id",
    "patient_group",
    "major_labl",
    "cell_type_original",
    "disease",
    "sex",
    "development_stage",
]

# Disease-progression ordinal for the trajectory evaluation.
PROGRESSION = {"myogenic": 0, "ischemic": 1, "fibrotic": 2}

print("Reading data")
adata = sc.read_h5ad(par["input"])
print(adata)

if adata.raw is None:
    raise ValueError("Expected raw counts in .raw for the CELLxGENE h5ad")

print("Rebuilding AnnData with raw counts in .X")
raw = adata.raw.to_adata()
obs = adata.obs[KEEP_OBS].copy()

# Numeric age from the CELLxGENE development_stage label, e.g.
# "44-year-old stage" -> 44.
obs["age"] = (
    obs["development_stage"].astype(str).str.extract(r"(\d+)").astype("float")
)

# Ordinal disease progression (myogenic -> ischemic -> fibrotic).
obs["disease_progression"] = obs["patient_group"].map(PROGRESSION).astype("Int64")

clean = ad.AnnData(X=raw.X, obs=obs, var=raw.var)

# Gene symbols as var_names; keep the Ensembl IDs in a column. Drop the
# now-redundant feature_name column: var_names_make_unique() can make the
# index differ from it, which anndata refuses to write.
if "feature_name" in clean.var.columns:
    clean.var["ensembl_id"] = clean.var_names
    clean.var_names = clean.var["feature_name"].astype(str)
    clean.var = clean.var.drop(columns="feature_name")
    clean.var_names_make_unique()
    clean.var.index.name = "gene_symbol"

print("Permissive QC filter")
sc.pp.filter_cells(clean, min_genes=200)
sc.pp.filter_genes(clean, min_cells=10)

# Cast to a plain CSR float32 of integer counts (drop float quirks).
clean.X = clean.X.astype(np.float32)

print(clean)
print("Samples:", clean.obs["sample"].nunique(), "| donors:", clean.obs["donor_id"].nunique())
print(clean.obs["patient_group"].value_counts(dropna=False))

print(f"Writing {par['output']} (compression={par.get('output_compression')})")
clean.write(par["output"], compression=par.get("output_compression"))
