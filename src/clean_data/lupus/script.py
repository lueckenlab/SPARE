import scanpy as sc
import anndata as ad
import requests
import pandas as pd

## VIASH START
par = {
    "input": "data/lupus.h5ad",
    "output": "data/lupus_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

print("Downloading clinical variables")
clinical_variables_url = "https://ucsf.app.box.com/shared/static/vbdcqo28ypy7wvgfoxti3blujub6zsaz"
response = requests.get(clinical_variables_url)
if response.status_code != 200:
    raise Exception(f"Failed to download clinical variables: {response.status_code}, {response.text}")
clinical_variables_df = pd.read_excel(response.content).drop(columns=["age", "female", "agedx", "raceeth", "raceethdesc", "subjectid", "X"])

# Load input adata
input_adata = sc.read_h5ad(par["input"])

# Merge clinical variables
obs = pd.merge(input_adata.obs, clinical_variables_df, left_on="ind_cov", right_on="genotypeid", how="left")

# Final adata
adata = ad.AnnData(var=input_adata.raw.var, obs=obs, layers={"X_raw_counts": input_adata.raw.X})

print("Writing output")
print(adata)
adata.write(par["output"], compression=par["output_compression"])