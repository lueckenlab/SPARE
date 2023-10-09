import scanpy as sc

adata = sc.read_h5ad("../data/HLCA_subset.h5ad")

print("Calculating HVG")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, flavor="seurat_v3", n_top_genes=3000)
adata = adata[:, adata.var.highly_variable].copy()

sc.tl.pca(adata)

adata.write("../data/HLCA_processed.h5ad")
