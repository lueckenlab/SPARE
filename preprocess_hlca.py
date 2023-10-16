from patient_representation.pp import is_count_data
import scanpy as sc

adata = sc.read_h5ad("../data/HLCA_subset.h5ad")

SAMPLE_KEY = "donor_id"
CELL_TYPE_KEY = "cell_type"
BATCH_KEY = "dataset"

print("X contains count data:", is_count_data(adata.X))
print("raw.X contains count data:", is_count_data(adata.raw.X))

# Copy raw counts to obsm
adata.obsm["raw"] = adata.raw.X.copy()
adata.layers["raw"] = adata.raw.X.copy()

print("Calculating HVG")
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000, batch_key=BATCH_KEY, layer="raw")

print("Normalizing data")
sc.pp.normalize_total(adata, target_sum=1e4)
print("Log-transforming data")
sc.pp.log1p(adata)

adata = adata[:, adata.var.highly_variable].copy()
adata.layers["raw"] = adata.layers["raw"][:, adata.index]
print("adata.shape", adata.shape)
print("adata.layers['raw'].shape", adata.layers["raw"].shape)
adata.obsm["raw"] = adata.layers["raw"]

sc.tl.pca(adata)

adata.write("../data/HLCA_processed.h5ad")
