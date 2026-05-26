import scanpy as sc
import anndata as ad
import requests
import pandas as pd


# Sample-level derived columns for evaluation. We expose them on the cell-level
# obs so they propagate through preprocess → metadata.csv → aggregate_representations
# without further wiring.
_DISEASE_STATUS_MAP = {
    "systemic lupus erythematosus": "SLE",
    "normal": "control",
}

_FLARE_STATUS_MAP = {
    # disease_state is already merged in by the upstream Perez et al. clinical
    # spreadsheet; we just rename for clarity.
    "flare": "pre-treatment",
    "treated": "post-treatment",
    # "managed" (stable SLE) and "na" (HC) intentionally drop to NaN — they're
    # not part of the flare/post-flare comparison.
}


def _derive_evaluation_columns(obs):
    """Add disease_status / flare_status / replicate / age columns in-place."""
    obs["disease_status"] = obs["disease"].astype(object).map(_DISEASE_STATUS_MAP)
    obs["flare_status"] = obs["disease_state"].astype(object).map(_FLARE_STATUS_MAP)

    # Replicate: True if this sample is not the first occurrence of its ind_cov.
    # First occurrence per donor is treated as the primary sample; downstream
    # KNN / trajectory metrics can filter out replicate=True to avoid double-
    # counting the same donor.
    sample_ids = obs["ind_cov_batch_cov"].astype(str)
    donor_ids = obs["ind_cov"].astype(str)
    sample_to_donor = (
        pd.DataFrame({"sample": sample_ids, "donor": donor_ids})
        .drop_duplicates("sample")
    )
    primary_sample_per_donor = sample_to_donor.drop_duplicates("donor")["sample"]
    primary = set(primary_sample_per_donor.tolist())
    obs["replicate"] = ~sample_ids.isin(primary)

    # Age: parse leading integer from "NN-year-old stage".
    obs["age"] = (
        obs["development_stage"]
        .astype(str)
        .str.extract(r"^(\d+)", expand=False)
        .astype(float)
    )
    return obs

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

# Recreate the Perez et al. sample identifier (donor × processing cohort).
# The current CXG redistribution drops this column, but the canonical sample
# definition counts the same donor processed in different cohorts as separate
# samples.
obs["ind_cov_batch_cov"] = obs["ind_cov"].astype(str) + "_" + obs["Processing_Cohort"].astype(str)

# pd.merge demotes categorical columns to object. Restore the cell-type
# column to category before downstream preprocess. Also fillna('nan')
# first - scvi-tools' SCANVI uses unlabeled_category='nan' and treats
# pandas NaN as -1 codes which trips its categorical-registration check
# even when the visible categories match.
obs["ct_cov"] = obs["ct_cov"].astype(object).fillna("nan").astype("category")
obs["ind_cov_batch_cov"] = obs["ind_cov_batch_cov"].astype("category")

obs = _derive_evaluation_columns(obs)

# Final adata
adata = ad.AnnData(var=input_adata.raw.var, obs=obs, layers={"X_raw_counts": input_adata.raw.X})

print("Writing output")
print(adata)
adata.write(par["output"], compression=par["output_compression"])