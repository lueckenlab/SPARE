import scanpy as sc

adata = sc.read_h5ad("../data/HLCA_subset.h5ad")

SAMPLE_KEY = "donor_id"
CELL_TYPE_KEY = "cell_type"
BATCH_KEY = "dataset"

# Copy raw counts to obsm
adata.obsm["raw"] = adata.X
adata.layers["raw"] = adata.X

print("Normalizing data")
sc.pp.normalize_total(adata, target_sum=1e4)
print("Log-transforming data")
sc.pp.log1p(adata)

print("Calculating HVG")
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000, batch_key=BATCH_KEY, layer="raw")

adata = adata[:, adata.var.highly_variable].copy()

sc.tl.pca(adata)

adata.write("../data/HLCA_processed.h5ad")
