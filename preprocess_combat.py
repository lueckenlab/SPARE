import scanpy as sc
import scvi
import pandas as pd

ADATA_PATH = "../data/COMBAT-CITESeq-DATA.h5ad"
IFN_1_SIGNATURE_PATH = "../data/ifn_1_score.tsv"

# Read IFN 1 signature genes
column_names = ["Gene/product", "gene_name", "Gene/product name", "Annotation qualifier",
                "GO class (direct)", "Annotation extension", "Contributor",
                "Organism", "Evidence", "Evidence with", "PANTHER family",
                "Type", "Isoform", "Reference", "Date"]
ifn_1_signature = pd.read_csv(IFN_1_SIGNATURE_PATH, sep="\t", names=column_names)
ifn_1_signature_genes = ifn_1_signature["gene_name"].unique()
print("Prepared ", len(ifn_1_signature_genes), "genes for calculating IFN1 signature")

print("Reading data")
adata = sc.read_h5ad(ADATA_PATH)

# COMBAT data is multimodal i.e. it contains protein expression as well
# We will leave it behind for now and focus only on the RNA expression data
print("Subsetting gene expression data")
adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()

# Remove cells with no label
adata = adata[adata.obs["Annotation_major_subset"] != "nan"]

# Find highly-variable genes
print("Subsetting HVG")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, flavor="seurat_v3",
                            n_top_genes=3000, layer="raw")

adata = adata[:, adata.var.highly_variable].copy()

print("Calculating IFN1 signature")
sc.tl.score_genes(adata, ifn_1_signature_genes, score_name="ifn_1_score")

# Run PCA
print("Running PCA")
sc.tl.pca(adata)

# Obtain scVI embedding
print("Running scVI")
adata.raw = adata
scvi.model.SCVI.setup_anndata(adata, layer="raw", batch_key="scRNASeq_sample_ID")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()

# Obtain scanVI embedding
print("Running scANVI")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="Annotation_major_subset",
    unlabeled_category="nan",
)
lvae.train(max_epochs=20, n_samples_per_label=100)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

# Use Pool_ID as a batch now
print("Running scVI with Pool_ID as a batch")
scvi.model.SCVI.setup_anndata(adata, layer="raw", batch_key="Pool_ID")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
adata.obsm["X_scVI_Pool_ID"] = vae.get_latent_representation()

# Obtain scanVI embedding
print("Running scANVI with Pool_ID as a batch")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="Annotation_major_subset",
    unlabeled_category="nan",
)
lvae.train(max_epochs=20, n_samples_per_label=100)
adata.obsm["X_scANVI_Pool_ID"] = lvae.get_latent_representation(adata)

adata.write("../data/combat_processed.h5ad")
