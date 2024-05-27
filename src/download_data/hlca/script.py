import boto3
from tqdm import tqdm
from botocore.config import Config
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError, BotoCoreError
from botocore import UNSIGNED
import cellxgene_census
import urllib.parse
import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
from patient_representation.pp import is_count_data



## VIASH START
par = {
    "output": "data/hlca.h5ad",
    "output_compression": "gzip",
    # "max_chunks": None -> this was complicated to implement with boto3
}
## VIASH END

HLCA_CELLXGENE_ID = "9f222629-9e39-47d0-b83f-e08d610c7479"

try:
    locator = cellxgene_census.get_source_h5ad_uri(dataset_id=HLCA_CELLXGENE_ID)
    uri = locator['uri']
    parsed_uri = urllib.parse.urlparse(uri)
    bucket = parsed_uri.netloc
    key = parsed_uri.path.lstrip('/')
    region_name = locator.get('s3_region')

    print("+++++++++++")
    print("URI:", uri)
    print("Bucket:", bucket)
    print("Key:", key)
    print("Region:", region_name)
    print("+++++++++++")
except Exception as e:
    print("Error retrieving file locator information:", e)
    sys.exit(1)

def download_file_from_s3(bucket, key, download_path, region_name):
    try:
        s3 = boto3.client('s3', region_name=region_name, config=Config(signature_version=UNSIGNED))
        response = s3.head_object(Bucket=bucket, Key=key)
        total_size_in_bytes = response['ContentLength']

        with tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True, desc='Downloading') as progress:
            def hook(bytes_transferred):
                progress.update(bytes_transferred)

            with open(download_path, 'wb') as f:
                s3.download_fileobj(bucket, key, f, Callback=hook)
        print("Download completed successfully.")

    except FileNotFoundError:
        print("Error: The file path provided does not exist.")
    except NoCredentialsError:
        print("Error: AWS credentials not available.")
    except PartialCredentialsError:
        print("Error: Incomplete AWS credentials.")
    except ClientError as e:
        print(f"Error: AWS client error - {e}")
    except BotoCoreError as e:
        print(f"Error: Low-level boto3 error - {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


download_dir = Path(par["output"]).parent

if not download_dir.exists():
    download_dir.mkdir(parents=True)

download_file_from_s3(bucket, key, par["output"], region_name)

# The whole list of diseases in HLCA. The ones that are not pneumonia-related are commented out.
pneumonia_related_diseases = [
    "disease",
    "normal",
    "COVID-19",
    "pulmonary fibrosis"
    "interstitial lung disease",
    # chronic obstructive pulmonary disease
    # lung adenocarcinoma
    "pneumonia",
    # chronic rhinitis
    # lung large cell carcinoma
    # squamous cell lung carcinoma
    # cystic fibrosis
    # lymphangioleiomyomatosis
    # pleomorphic carcinoma
    # hypersensitivity pneumonitis
    "non-specific interstitial pneumonia",
    # pulmonary sarcoidosis
]

adata = sc.read_h5ad(par["output"])

print("HLCA shape before filtering out diseases:", adata.shape)
adata = adata[adata.obs["disease"].isin(pneumonia_related_diseases)]
print("HLCA shape after filtering out diseases:", adata.shape)

print("Filtering cells with no annotation")
adata = adata[pd.notna(adata.obs["study"])]
print("HLCA shape after filtering out unannotated cells:", adata.shape)

print("X contains count data:", is_count_data(adata.X))
print("raw.X contains count data:", is_count_data(adata.raw.X))

# Copy raw counts to obsm
#but whas is need of Copy() and having in both layers and obsm
#maybe at least delete raw.X after first copy and then start next copy step
adata.obsm["raw"] = adata.raw.X.copy()
adata.layers["raw"] = adata.raw.X.copy()

adata.write(par["output"], compression=par["output_compression"])
