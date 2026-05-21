"""Clean the SEA-AD Whole Taxonomy datasets (MTG and DLPFC share this).

SEA-AD = Seattle Alzheimer's Disease Brain Cell Atlas (Gabitto et al.
2024). The CELLxGENE files (schema 7.0.0) carry normalized counts in
`.X`, raw integer counts in `.raw`, full-gene Ensembl var with symbols in
`var["feature_name"]`. ~1.4M nuclei x 35k genes per region.

Steps (memory-conscious — these files are ~50GB on disk):
1. Pull raw counts out of `.raw` into `.X`; drop the normalized matrix
   immediately so we never hold both copies after the read.
2. Gene symbols as var_names (Ensembl kept in a column).
3. Curate a donor-level-meaningful obs and derive ordinal pathology
   scores for the AD continuum (adnc_score is the trajectory variable;
   braak/cerad/thal also provided). Neurotypical "Reference" donors map
   to 0 (no AD neuropathologic change).
4. Permissive QC.
"""
import gc
import re

import anndata as ad
import numpy as np
import scanpy as sc

## VIASH START
par = {
    "input": "data/sea_ad_mtg/sea_ad_mtg.h5ad",
    "output": "data/sea_ad_mtg/sea_ad_mtg_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

KEEP_OBS = [
    "donor_id",
    "Specimen ID",
    "Class",
    "Subclass",
    "Supertype",
    "assay",
    "disease",
    "sex",
    "Neurotypical reference",
    "Cognitive status",
    "ADNC",
    "Braak stage",
    "Thal phase",
    "CERAD score",
    "APOE4 status",
    "Age at death",
    "Years of education",
    "PMI",
    "Lewy body disease pathology",
    "LATE-NC stage",
    "Microinfarct pathology",
    "self_reported_ethnicity",
    "development_stage",
]

# Ordinal encodings of the AD-continuum neuropathology. "Reference"
# (neurotypical reference donors) -> 0, i.e. no AD neuropathologic change.
ADNC_MAP = {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3, "Reference": 0}
BRAAK_MAP = {
    "Braak 0": 0, "Braak I": 1, "Braak II": 2, "Braak III": 3,
    "Braak IV": 4, "Braak V": 5, "Braak VI": 6, "Reference": 0,
}
CERAD_MAP = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3, "Reference": 0}
THAL_MAP = {f"Thal {i}": i for i in range(6)}
THAL_MAP["Reference"] = 0


def ordinal(series, mapping, name):
    mapped = series.map(mapping)
    unmapped = sorted(set(series.dropna().unique()) - set(mapping))
    if unmapped:
        print(f"  WARNING: unmapped {name} levels (-> NaN): {unmapped}")
    return mapped.astype("Int64")


print("Reading data")
adata = sc.read_h5ad(par["input"])
print(adata)

if adata.raw is None:
    raise ValueError("Expected raw counts in .raw for the CELLxGENE h5ad")

print("Extracting raw counts into .X; dropping normalized matrix")
raw = adata.raw.to_adata()                       # raw counts, full gene set
obs = adata.obs[KEEP_OBS].copy()
del adata
gc.collect()

# Derived ordinal trajectory scores.
obs["adnc_score"] = ordinal(obs["ADNC"], ADNC_MAP, "ADNC")
obs["braak_score"] = ordinal(obs["Braak stage"], BRAAK_MAP, "Braak stage")
obs["cerad_score"] = ordinal(obs["CERAD score"], CERAD_MAP, "CERAD score")
obs["thal_score"] = ordinal(obs["Thal phase"], THAL_MAP, "Thal phase")

# Numeric age from development_stage ("84-year-old stage" -> 84,
# "80 year-old and over stage" -> 80).
obs["age"] = obs["development_stage"].astype(str).str.extract(r"(\d+)").astype("float")

clean = ad.AnnData(X=raw.X, obs=obs, var=raw.var)
del raw
gc.collect()

# Gene symbols as var_names; keep Ensembl IDs in a column. Drop the
# feature_name column so the made-unique index can't collide with it.
if "feature_name" in clean.var.columns:
    clean.var["ensembl_id"] = clean.var_names
    clean.var_names = clean.var["feature_name"].astype(str)
    clean.var = clean.var.drop(columns="feature_name")
    clean.var_names_make_unique()
    clean.var.index.name = "gene_symbol"

print("Permissive QC filter")
sc.pp.filter_cells(clean, min_genes=200)
sc.pp.filter_genes(clean, min_cells=10)
clean.X = clean.X.astype(np.float32)

print(clean)
print("Samples (donors):", clean.obs["donor_id"].nunique(),
      "| Subclasses:", clean.obs["Subclass"].nunique())
print("ADNC:\n", clean.obs.groupby("donor_id", observed=True)["adnc_score"].first().value_counts(dropna=False).sort_index())

print(f"Writing {par['output']} (compression={par.get('output_compression')})")
clean.write(par["output"], compression=par.get("output_compression"))
