"""Split the combined cleaned cohort into two incomparable datasets by CD45
sort: CD45+ (immune) and CD45- (tumor/stroma). Within each split `sort` is
constant and the source-dataset split encodes the response label, so there is
no batch covariate; integration will be over sample_id only.
"""
import os, scanpy as sc

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
src = os.path.join(HERE, "scct_db_cleaned.h5ad")
adata = sc.read_h5ad(src)
print("combined:", adata.shape, "| sort:", adata.obs["sort"].value_counts().to_dict())

# Drop unannotated cells (the only junk label in this cohort is "Not assigned").
n0 = adata.n_obs
adata = adata[adata.obs["cell_type"] != "Not assigned"].copy()
print(f"dropped {n0 - adata.n_obs} 'Not assigned' cells -> {adata.n_obs} cells")

SPLITS = {
    "scct_db_cd45pos": "CD45+",
    "scct_db_cd45neg": "CD45-",
}
for ds_id, sort_val in SPLITS.items():
    sub = adata[adata.obs["sort"] == sort_val].copy()
    # drop genes now empty in this compartment
    sc.pp.filter_genes(sub, min_cells=10)
    out_dir = os.path.join(REPO, "data", ds_id)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{ds_id}_cleaned.h5ad")
    sub.write(out, compression="gzip")
    print(f"{ds_id}: {sub.shape} | samples={sub.obs['sample_id'].nunique()} "
          f"| response={sub.obs['response'].value_counts().to_dict()} -> {out}")
