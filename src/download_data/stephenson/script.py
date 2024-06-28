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

print("Copying raw counts to layers")
adata.X = adata.raw.X.copy()

try:
    print("X contains count data:", is_count_data(adata.X))
except Exception as e:
    print("Error checking if adata.X contains count data:", e)

decade_mapping = {
	'sixth decade human stage': '50',
	'fourth decade human stage': '30',
	'fifth decade human stage': '40',
	'third decade human stage': '20',
	'eighth decade human stage': '70',
	'seventh decade human stage': '60',
	'ninth decade human stage': '80',
	'49-year-old human stage': '40',  # Assuming this falls into the fourth decade
	'52-year-old human stage': '50',  # Assuming this falls into the fifth decade
	'70-year-old human stage': '70',  # Assuming this falls into the seventh decade
	'66-year-old human stage': '60',  # Assuming this falls into the sixth decade
	'58-year-old human stage': '50',
	'63-year-old human stage': '60',
	'51-year-old human stage': '50',
	'65-year-old human stage': '60',
	'26-year-old human stage': '20',
	'30-year-old human stage': '30',
	'38-year-old human stage': '30',
	'77-year-old human stage': '70',
	'76-year-old human stage': '70',
	'73-year-old human stage': '70',
	'71-year-old human stage': '70',
	'50-year-old human stage': '50',
	'69-year-old human stage': '60',
	'44-year-old human stage': '40',
	'64-year-old human stage': '60',
	'54-year-old human stage': '50',
	'55-year-old human stage': '50',
	'62-year-old human stage': '60',
	'56-year-old human stage': '50',
	'47-year-old human stage': '40',
	'57-year-old human stage': '50',
	'87-year-old human stage': '80',
	'31-year-old human stage': '30',
	'39-year-old human stage': '30',
	'40-year-old human stage': '40',
	'92-year-old human stage': '90',
	'83-year-old human stage': '80',
	'84-year-old human stage': '80',
	'46-year-old human stage': '40',
	'80-year-old human stage': '80',
	'59-year-old human stage': '50',
	'60-year-old human stage': '60',
	'25-year-old human stage': '20',
	'68-year-old human stage': '60',
	'21-year-old human stage': '20',
}
# Map the values in the "development_stage" column to the corresponding decades
adata.obs['development_stage'] = adata.obs['development_stage'].map(decade_mapping)
adata.obs['development_stage'] = adata.obs['development_stage'].astype('int')

#Change the type of Collection day and merge the data since 4 classes with only one patient
mapping = {
    'D0': 0,
    'D28': 28,
    'D7': 10,
    'D9': 10,
    'D12': 10,
    'D13': 10,
}

adata.obs['Collection_Day'] = adata.obs['Collection_Day'].replace(mapping).astype(int)


adata.obs.loc[
    adata.obs["Days_from_onset"].isin(["Healthy", "Not_known", "LPS", "nan", "Non_covid"]),
    "Days_from_onset"] = np.nan

adata.obs["Days_from_onset"] = adata.obs["Days_from_onset"].cat.codes.astype(float)
# Properly set NaN values
adata.obs.loc[
    adata.obs["Days_from_onset"] == -1,
    "Days_from_onset"
] = np.nan

adata = adata[adata.obs['Status'] != 'LPS']
# adata = adata[adata.obs['Site'] != 'Sanger'] ## Sanger has 74913 cells

del adata.raw

adata.layers["X_raw_counts"] = adata.X.copy()
print("adata.obsm['X_raw_counts'].shape", adata.obsm["X_raw_counts"].shape)

#TODO: rest of the dataset specefic preprocessing
print(adata)
print("Saving output")
adata.write(par["output"], compression=par["output_compression"])