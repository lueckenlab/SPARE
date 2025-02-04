import pandas as pd
import scanpy as sc
from pathlib import Path

## VIASH START
par = {
    "input": "data/combat.h5ad",
    "output": "data/combat_cleaned.h5ad",
    "output_compression": "gzip",
}
## VIASH END

IFN_1_SIGNATURE_PATH = "../../../data/ifn_1_score.tsv"

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
adata = adata[adata.obs["Annotation_major_subset"] != "nan"].copy()

print("Calculating IFN1 signature")
sc.tl.score_genes(adata, ifn_1_signature_genes, score_name="ifn_1_score")

print(adata)

# Extract raw counts that are needed by some methods
print("adata.layers['raw'].shape", adata.layers["raw"].shape)
print("adata.shape", adata.shape)

print("Moving raw counts to X")
adata.X = adata.layers["raw"].copy()
del adata.layers["raw"]

print("Copying raw counts to layers")
adata.layers["X_raw_counts"] = adata.X.copy()

print("Saving output")
adata.write(par["output"], compression=par["output_compression"])
