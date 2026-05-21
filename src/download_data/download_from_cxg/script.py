"""Download a dataset from CELLxGENE Discover by its dataset ID.

Resolves the dataset's H5AD asset URL through the public Discover
curation API and streams it over HTTPS. This deliberately avoids
`cellxgene_census.get_source_h5ad_uri`: the census build pinned in this
environment rejects the current census release with
`ValueError: Unsupported SOMA object encoding version 1.1.0`. The asset
URLs the API returns (datasets.cellxgene.cziscience.com/...) are public
and need neither the census nor signed S3 access.
"""
import sys
from pathlib import Path

import requests
from tqdm import tqdm

## VIASH START
par = {
    "input": "1c739a3e-c3f5-49d5-98e0-73975e751201",
    "output": "output.h5ad",
    "max_chunks": None,
}
## VIASH END

CHUNK_SIZE = 1024 * 1024  # 1 MB
DATASETS_INDEX = "https://api.cellxgene.cziscience.com/curation/v1/datasets"


def resolve_h5ad_url(dataset_id):
    """Return the H5AD download URL for a CELLxGENE Discover dataset ID.

    The per-dataset endpoint requires the collection ID, which we do not
    have here, so we look the dataset up in the global datasets index;
    each entry carries its asset download URLs directly.
    """
    print(f"Resolving H5AD asset for dataset {dataset_id}")
    resp = requests.get(DATASETS_INDEX, timeout=120)
    resp.raise_for_status()
    datasets = resp.json()

    match = next((d for d in datasets if d.get("dataset_id") == dataset_id), None)
    if match is None:
        raise ValueError(
            f"Dataset ID {dataset_id!r} not found in the CELLxGENE Discover index "
            f"({len(datasets)} datasets searched)."
        )

    h5ad = next(
        (a for a in match.get("assets", []) if a.get("filetype") == "H5AD"), None
    )
    if h5ad is None or not h5ad.get("url"):
        raise ValueError(f"No H5AD asset with a download URL for dataset {dataset_id!r}.")
    return h5ad["url"]


def download(url, output, max_chunks=None):
    """Stream `url` to `output`, optionally stopping after max_chunks (1 MB each)."""
    Path(output).parent.mkdir(parents=True, exist_ok=True)

    with requests.get(url, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("Content-Length", 0)) or None
        if max_chunks is not None:
            total = min(total, max_chunks * CHUNK_SIZE) if total else max_chunks * CHUNK_SIZE

        with tqdm(total=total, unit="iB", unit_scale=True, desc="Downloading") as bar, \
                open(output, "wb") as f:
            for i, chunk in enumerate(resp.iter_content(chunk_size=CHUNK_SIZE)):
                f.write(chunk)
                bar.update(len(chunk))
                if max_chunks is not None and i + 1 >= max_chunks:
                    print(f"Stopped after {max_chunks} chunk(s) (max_chunks set).")
                    break
    print("Download completed successfully.")


try:
    url = resolve_h5ad_url(par["input"])
    print(f"H5AD URL: {url}")
    download(url, par["output"], max_chunks=par.get("max_chunks"))
except requests.HTTPError as e:
    print(f"Error: HTTP request to CELLxGENE failed - {e}")
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    sys.exit(1)
