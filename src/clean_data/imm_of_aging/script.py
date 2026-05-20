"""Clean the AIFI Immunobiology of Aging concat h5ad for SPARE.

Three things only:
1. Cast X from uint16 to int32 (scipy ops in src/preprocess reject uint16).
2. Parse `subject.ageAtFirstDraw` "89+" -> 89 so the column is numeric.
3. Drop per-cell high-cardinality categoricals that bloat gzip writes.

Then apply a permissive QC pass. Raw counts stay in `.X`; the
`X_raw_counts` layer is created in src/preprocess after HVG subsetting
where the matrix is 6x smaller (a full-resolution copy here would blow
past 400 GB on this cohort).
"""
import gc

import numpy as np
import pandas as pd
import scanpy as sc

## VIASH START
par = {
    "input": "data/imm_of_aging/imm_of_aging.h5ad",
    "output": "data/imm_of_aging/imm_of_aging_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

DROP_COLS = ["cell_name", "original_barcodes"]
AGE_COL = "subject.ageAtFirstDraw"

print("Reading data")
adata = sc.read_h5ad(par["input"])
print(adata)

if adata.X.dtype == np.uint16:
    print("Casting X uint16 -> int32 (in place)")
    adata.X.data = adata.X.data.astype(np.int32)
    gc.collect()

for col in DROP_COLS:
    if col in adata.obs.columns:
        print(f"Dropping high-cardinality column: {col}")
        adata.obs = adata.obs.drop(columns=col)

if AGE_COL in adata.obs.columns:
    raw = adata.obs[AGE_COL].astype(str).str.rstrip("+")
    adata.obs[AGE_COL] = pd.to_numeric(raw, errors="coerce")
    n_bad = int(adata.obs[AGE_COL].isna().sum())
    print(f"Parsed {AGE_COL} to numeric "
          f"(range [{adata.obs[AGE_COL].min()}, {adata.obs[AGE_COL].max()}], "
          f"{n_bad:,} unparseable cells)")

print("Permissive QC filter")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
print(adata)

print(f"Writing {par['output']} (compression={par.get('output_compression')})")
adata.write(par["output"], compression=par.get("output_compression"))
