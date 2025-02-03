import boto3
from tqdm import tqdm
from botocore.config import Config
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError, BotoCoreError
from botocore import UNSIGNED
import cellxgene_census
import urllib.parse
import sys
from pathlib import Path

## VIASH START
par = {
    "input": "3faad104-2ab8-4434-816d-474d8d2641db",
    "output": "output.h5ad",
    "max_chunks": None
}
## VIASH END

def download_file_from_s3(bucket, key, download_path, region_name, max_chunks=None, chunk_size=1024*1024):
    """Download a file from S3 with progress tracking.

    Parameters
    ----------
    bucket : str
        S3 bucket name
    key : str
        S3 object key
    download_path : str
        Path where to save the downloaded file
    region_name : str
        AWS region name
    max_chunks : int, optional
        Maximum number of chunks to download, by default None
    chunk_size : int, optional
        Size of chunks to download in bytes, by default 1MB

    Returns
    -------
    None
    """
    try:
        s3 = boto3.client("s3", region_name=region_name, config=Config(signature_version=UNSIGNED))
        response = s3.head_object(Bucket=bucket, Key=key)
        total_size_in_bytes = response["ContentLength"]

        with tqdm(total=total_size_in_bytes, unit="iB", unit_scale=True, desc="Downloading") as progress:
            def hook(bytes_transferred):
                progress.update(bytes_transferred)

            if max_chunks is not None:
                # If max_chunks is set, only download a portion of the file
                range_end = min(max_chunks * chunk_size, total_size_in_bytes)
                response = s3.get_object(
                    Bucket=bucket,
                    Key=key,
                    Range=f"bytes=0-{range_end-1}"
                )
                with open(download_path, "wb") as f:
                    for chunk in response["Body"].iter_chunks(chunk_size=chunk_size):
                        f.write(chunk)
                        progress.update(len(chunk))
            else:
                # Download the entire file
                with open(download_path, "wb") as f:
                    s3.download_fileobj(bucket, key, f, Callback=hook)

        print("Download completed successfully.")

    except FileNotFoundError:
        print("Error: The file path provided does not exist.")
        sys.exit(1)
    except NoCredentialsError:
        print("Error: AWS credentials not available.")
        sys.exit(1)
    except PartialCredentialsError:
        print("Error: Incomplete AWS credentials.")
        sys.exit(1)
    except ClientError as e:
        print(f"Error: AWS client error - {e}")
        sys.exit(1)
    except BotoCoreError as e:
        print(f"Error: Low-level boto3 error - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

try:
    locator = cellxgene_census.get_source_h5ad_uri(dataset_id=par["input"])
    uri = locator["uri"]
    parsed_uri = urllib.parse.urlparse(uri)
    bucket = parsed_uri.netloc
    key = parsed_uri.path.lstrip("/")
    region_name = locator.get("s3_region")

    print("+++++++++++")
    print("URI:", uri)
    print("Bucket:", bucket)
    print("Key:", key)
    print("Region:", region_name)
    print("+++++++++++")
except Exception as e:
    print("Error retrieving file locator information:", e)
    sys.exit(1)

download_dir = Path(par["output"]).parent

if not download_dir.exists():
    download_dir.mkdir(parents=True)

download_file_from_s3(bucket, key, par["output"], region_name, max_chunks=par["max_chunks"]) 
