import pandas as pd
import scanpy as sc
import numpy as np
from patient_representation.pp import is_count_data

## VIASH START
par = {
    "input": "data/stephenson.h5ad",
    "output": "data/stephenson_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

print("Reading data")
adata = sc.read_h5ad(par["input"])

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
    "sixth decade stage": "50",
    "fourth decade stage": "30",
    "fifth decade stage": "40",
    "third decade stage": "20",
    "eighth decade stage": "70",
    "seventh decade stage": "60",
    "ninth decade stage": "80",
    "49-year-old stage": "40",  # Assuming this falls into the fourth decade
    "52-year-old stage": "50",  # Assuming this falls into the fifth decade
    "70-year-old stage": "70",  # Assuming this falls into the seventh decade
    "66-year-old stage": "60",  # Assuming this falls into the sixth decade
    "58-year-old stage": "50",
    "63-year-old stage": "60",
    "51-year-old stage": "50",
    "65-year-old stage": "60",
    "26-year-old stage": "20",
    "30-year-old stage": "30",
    "38-year-old stage": "30",
    "77-year-old stage": "70",
    "76-year-old stage": "70",
    "73-year-old stage": "70",
    "71-year-old stage": "70",
    "50-year-old stage": "50",
    "69-year-old stage": "60",
    "44-year-old stage": "40",
    "64-year-old stage": "60",
    "54-year-old stage": "50",
    "55-year-old stage": "50",
    "62-year-old stage": "60",
    "56-year-old stage": "50",
    "47-year-old stage": "40",
    "57-year-old stage": "50",
    "87-year-old stage": "80",
    "31-year-old stage": "30",
    "39-year-old stage": "30",
    "40-year-old stage": "40",
    "92-year-old stage": "90",
    "83-year-old stage": "80",
    "84-year-old stage": "80",
    "46-year-old stage": "40",
    "80-year-old stage": "80",
    "59-year-old stage": "50",
    "60-year-old stage": "60",
    "25-year-old stage": "20",
    "68-year-old stage": "60",
    "21-year-old stage": "20",
}
# Map the values in the "development_stage" column to the corresponding decades
adata.obs["development_stage"] = adata.obs["development_stage"].map(decade_mapping)
adata.obs["development_stage"] = adata.obs["development_stage"].astype("int")

#Change the type of Collection day and merge the data since 4 classes with only one patient
mapping = {
    "D0": 0,
    "D28": 28,
    "D7": 10,
    "D9": 10,
    "D12": 10,
    "D13": 10,
}

adata.obs["Collection_Day"] = adata.obs["Collection_Day"].replace(mapping).astype(int)


adata.obs.loc[
    adata.obs["Days_from_onset"].isin(["Healthy", "Not_known", "LPS", "nan", "Non_covid"]),
    "Days_from_onset"] = np.nan

adata.obs["Days_from_onset"] = adata.obs["Days_from_onset"].cat.codes.astype(float)
# Properly set NaN values
adata.obs.loc[
    adata.obs["Days_from_onset"] == -1,
    "Days_from_onset"
] = np.nan

adata = adata[adata.obs["Status"] != "LPS"].copy()
# adata = adata[adata.obs["Site"] != "Sanger"] ## Sanger has 74913 cells

# Map severity to numeric values to calculate correlations
adata.obs["Status_on_day_collection_summary"] = adata.obs["Status_on_day_collection_summary"].map(
    {
        "Severe": 6.0,
        "Critical": 5.0,
        "Moderate": 4.0,
        "Mild": 3.0,
        "Asymptomatic": 2.0,
        "Healthy": 1.0,
        "Non_covid": np.nan
    }
)

del adata.raw

adata.layers["X_raw_counts"] = adata.X.copy()
print("adata.obsm['X_raw_counts'].shape", adata.layers["X_raw_counts"].shape)

#TODO: rest of the dataset specefic preprocessing
print(adata)
print("Saving output")
adata.write(par["output"], compression=par["output_compression"])