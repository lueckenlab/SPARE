import requests
from pathlib import Path

import pandas as pd
import scanpy as sc
from tqdm import tqdm
import numpy as np
from patient_representation.pp import is_count_data

## VIASH START
par = {
    "output": "data/stephenson.h5ad",
    "output_compression": "gzip",
    "max_chunks": None
}
## VIASH END

STEPHENSON_URL = "https://datasets.cellxgene.cziscience.com/9a305865-e7b0-4e8a-bbbf-0b63287e85fd.h5ad"

def download_file(url, filename, max_chunks=None, block_size = 131072):
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        total_size_in_bytes = int(response.headers.get('content-length', 0))
    
        with open(filename, "wb") as file:
            with tqdm(unit="iB", unit_scale=True, total=total_size_in_bytes) as progress:
                for chunk_number, data in enumerate(response.iter_content(block_size)):
                    if max_chunks is not None and chunk_number >= max_chunks:
                        break

                    file.write(data)
                    progress.update(len(data))


download_dir = Path(par["output"]).parent

if not download_dir.exists():
    download_dir.mkdir(parents=True)

download_file(STEPHENSON_URL, par["output"], max_chunks=par["max_chunks"])

print("Reading data")
adata = sc.read_h5ad(par["output"])

try:
    print("X contains count data:", is_count_data(adata.X))
except Exception as e:
    print("Error checking if adata.X contains count data:", e)

try:
    print("raw.X contains count data:", is_count_data(adata.raw.X))
except Exception as e:
    print("Error checking if raw.X contains count data:", e)

adata.obs.loc[
    adata.obs["Days_from_onset"].isin(["Healthy", "Not_known", "LPS", "nan", "Non_covid"]),
    "Days_from_onset"] = np.nan

adata.obs["Days_from_onset"] = adata.obs["Days_from_onset"].cat.codes.astype(float)
# Properly set NaN values
adata.obs.loc[
    adata.obs["Days_from_onset"] == -1,
    "Days_from_onset"
] = np.nan


print("Copying raw counts to layers")
adata.X = adata.raw.X.copy()
del adata.raw
adata.obsm["X_raw_counts"]  = adata.X.copy()
adata.layers["X_raw_counts"] = adata.X.copy()
print("adata.obsm['X_raw_counts'].shape", adata.obsm["X_raw_counts"].shape)

print(adata)
print("Saving output")
adata.write(par["output"], compression=par["output_compression"])