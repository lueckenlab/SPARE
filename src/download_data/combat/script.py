import requests
from pathlib import Path

import pandas as pd
import scanpy as sc
from tqdm import tqdm


## VIASH START
par = {
    "output": "combat.h5mu",
    "output_compression": "gzip",
    "max_chunks": None
}
## VIASH END

COMBAT_URL = "https://zenodo.org/record/6120249/files/COMBAT-CITESeq-DATA.h5ad"
IFN_1_SIGNATURE_PATH = "data/ifn_1_score.tsv"

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
adata = adata[adata.obs["Annotation_major_subset"] != "nan"]

adata.write(par["output"], compression=par["output_compression"])
