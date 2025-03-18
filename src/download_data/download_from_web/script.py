import requests
from pathlib import Path
from tqdm import tqdm

## VIASH START
par = {
    "input": "https://example.com/file.h5ad",
    "output": "output.h5ad",
    "max_chunks": None
}
## VIASH END

def download_file(url, filename, max_chunks=None, block_size=131072):
    """Download a file from a URL with progress bar.

    Parameters
    ----------
    url : str
        URL to download from
    filename : str
        Path where to save the downloaded file
    max_chunks : int, optional
        Maximum number of chunks to download, by default None
    block_size : int, optional
        Size of chunks to download in bytes, by default 131072

    Returns
    -------
    None
    """
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        total_size_in_bytes = int(response.headers.get("content-length", 0))
    
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

download_file(par["input"], par["output"], max_chunks=par["max_chunks"]) 