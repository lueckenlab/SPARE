import requests
from pathlib import Path

from tqdm import tqdm


## VIASH START
par = {
    "output": "combat.h5mu",
    "output_compression": "gzip",
    "max_chunks": None
}
## VIASH END

COMBAT_URL = "https://zenodo.org/record/6120249/files/COMBAT-CITESeq-DATA.h5ad"

def download_file(url, filename, max_chunks=None):
    response = requests.get(url, stream=True)

    with open(filename, "wb") as f:
        for chunk_number, data_chunk in enumerate(tqdm(response.iter_content())):
            f.write(data_chunk)

            if max_chunks is not None and chunk_number >= max_chunks:
                break


download_dir = Path(par["output"]).parent

if not download_dir.exists():
    download_dir.mkdir(parents=True)

download_file(COMBAT_URL, par["output"], max_chunks=par["max_chunks"])
