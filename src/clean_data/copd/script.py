import scanpy as sc

## VIASH START
par = {
    "input": "data/copd.h5ad",
    "output": "data/copd_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

print("Reading data")
adata = sc.read_h5ad(par["input"])

# Map GOLD_Stage to numeric values to calculate correlations
adata.obs["GOLD_Stage"] = adata.obs["GOLD_Stage"].map(
    {
        "Control": 0.0,
        "GOLD0": 1.0,
        "GOLD1_2": 2.0,
        "GOLD3_4": 3.0,
    }
)

print(adata)
print("Saving output")
adata.write(par["output"], compression=par["output_compression"])