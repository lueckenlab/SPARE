import scanpy as sc
import scvi

ADATA_PATH = "../data/onek1k.h5ad"

print("Reading data")
adata = sc.read_h5ad(ADATA_PATH)

SAMPLE_KEY = "donor_id"
CELL_TYPE_KEY = "cell_type"
BATCH_KEY = "pool_number"

# Copy raw counts to obsm
adata.obsm["raw"] = adata.X

print("Normalizing data")
sc.pp.normalize_total(adata, target_sum=1e4)
print("Log-transforming data")
sc.pp.log1p(adata)

# Find highly-variable genes
print("Subsetting HVG")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, flavor="seurat_v3",
                            n_top_genes=3000, layer="raw")
adata = adata[:, adata.var.highly_variable].copy()

# Obtain scVI embedding
print("Running scVI")
adata.raw = adata
scvi.model.SCVI.setup_anndata(adata, layer="raw", batch_key=BATCH_KEY)
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()

# Obtain scanVI embedding
print("Running scANVI")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key=CELL_TYPE_KEY,
    unlabeled_category="nan",
)
lvae.train(max_epochs=20, n_samples_per_label=100)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

adata.write("../data/onek1k_processed.h5ad")
