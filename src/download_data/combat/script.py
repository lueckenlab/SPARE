import requests
from pathlib import Path

import pandas as pd
import scanpy as sc
from tqdm import tqdm


## VIASH START
par = {
    "output": "combat.h5ad",
    "output_compression": "gzip",
    "max_chunks": None
}
## VIASH END

# TODO: check out downloading combat from cellxgene LATER
COMBAT_URL = "https://zenodo.org/record/6120249/files/COMBAT-CITESeq-DATA.h5ad"
IFN_1_SIGNATURE_PATH = "data/ifn_1_score.tsv"

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

download_file(COMBAT_URL, par["output"], max_chunks=par["max_chunks"])


# Read IFN 1 signature genes
column_names = ("Gene/product", "gene_name", "Gene/product name", "Annotation qualifier",
                "GO class (direct)", "Annotation extension", "Contributor",
                "Organism", "Evidence", "Evidence with", "PANTHER family",
                "Type", "Isoform", "Reference", "Date")

ifn_1_signature = pd.read_csv(IFN_1_SIGNATURE_PATH, sep="\t", names=column_names)
ifn_1_signature_genes = ifn_1_signature["gene_name"].unique()
print("Prepared ", len(ifn_1_signature_genes), "genes for calculating IFN1 signature")

print("Reading data")
adata = sc.read_h5ad(par["output"])

# COMBAT data is multimodal i.e. it contains protein expression as well
# We will leave it behind for now and focus only on the RNA expression data
print("Subsetting gene expression data")
is_rna_expression = adata.var["feature_types"] == "Gene Expression"
adata = adata[:, is_rna_expression].copy()

# Remove cells with no label
print("Removing cells with no label")
adata = adata[adata.obs["Annotation_major_subset"] != "nan"]

# Extract raw counts that are needed by some methods
print("Moving raw counts to obsm and layers")
adata.layers["raw"] = adata.layers["raw"][:, is_rna_expression]
print("adata.layers['raw'].shape", adata.layers["raw"].shape)

adata.X = adata.layers["raw"].copy()
del adata.layers["raw"]

print("Copying raw counts to layers")
adata.layers["X_raw_counts"] = adata.X.copy()

print("adata.obsm['X_raw_counts'].shape", adata.obsm["X_raw_counts"].shape)
print("Calculating IFN1 signature")
sc.tl.score_genes(adata, ifn_1_signature_genes, score_name="ifn_1_score")

print(adata)

print("Saving output")
adata.write(par["output"], compression=par["output_compression"])
