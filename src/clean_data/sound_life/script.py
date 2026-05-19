"""Clean the AIFI Sound Life concat h5ad for SPARE.

Two things only:
1. Cast X from uint16 to int32 (scipy ops in src/preprocess reject uint16).
2. Drop per-cell high-cardinality categoricals that bloat gzip writes.

Then stash raw counts in `layers["X_raw_counts"]` and apply a permissive
QC pass. The AIFI release ships pre-filtered; the filter is a safety net.
"""
import gc

import numpy as np
import scanpy as sc

## VIASH START
par = {
    "input": "data/sound_life/sound_life.h5ad",
    "output": "data/sound_life/sound_life_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

DROP_COLS = ["cell_name", "original_barcodes"]

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

print("Stashing raw counts in layers['X_raw_counts']")
adata.layers["X_raw_counts"] = adata.X.copy()

print("Permissive QC filter")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
print(adata)

print(f"Writing {par['output']} (compression={par.get('output_compression')})")
adata.write(par["output"], compression=par.get("output_compression"))
