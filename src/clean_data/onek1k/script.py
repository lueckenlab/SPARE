import scanpy as sc
from patient_representation.pp import is_count_data

## VIASH START
par = {
    "input": "data/onek1k.h5ad",
    "output": "data/onek1k_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

adata = sc.read_h5ad(par["input"])

print("X contains count data:", is_count_data(adata.X))

# Copy raw counts to obsm
adata.layers["X_raw_counts"] = adata.X.copy()

print(adata)
adata.write(par["output"], compression=par["output_compression"])