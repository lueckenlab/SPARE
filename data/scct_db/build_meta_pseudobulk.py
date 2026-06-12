"""Enrich the sample-level meta-AnnData (one obs = one sample; representations
already in .obsm) with:
  - obsm['cell_type_proportions']         (n_samples x n_cell_types)
  - obsm['pseudobulk_<cell_type>']        (n_samples x n_genes), one per cell
    type = mean log-norm expression over that sample's cells of that type
    (NaN where a sample has no cells of the type).
Gene order + cell-type maps are stored in .uns.

Usage: python build_meta_pseudobulk.py scct_db_cd45pos
"""
import os, sys, re, numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc

REPO = "/ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark"
ds = sys.argv[1] if len(sys.argv) > 1 else "scct_db_cd45pos"

meta = sc.read_h5ad(os.path.join(REPO, "output_scct_db", f"{ds}.sample_representation.output.h5ad"))
proc = sc.read_h5ad(os.path.join(REPO, "data", ds, f"{ds}_processed.h5ad"))
X = proc.layers["X_log_norm_counts"]
X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
genes = proc.var_names.tolist()
samp = proc.obs["sample_id"].astype(str).values
ct = proc.obs["cell_type"].astype(str).values

# --- per (sample, cell_type) mean of log-norm expression via sparse matmul ---
gkey = np.array([f"{s}|||{c}" for s, c in zip(samp, ct)])
uniq, inv = np.unique(gkey, return_inverse=True)
counts = np.bincount(inv, minlength=len(uniq)).astype(float)
G = sp.csr_matrix((1.0 / counts[inv], (inv, np.arange(X.shape[0]))), shape=(len(uniq), X.shape[0]))
PB = np.asarray((G @ X).todense())          # n_groups x n_genes
gsamp = np.array([u.split("|||")[0] for u in uniq])
gct = np.array([u.split("|||")[1] for u in uniq])
group_of = {(s, c): i for i, (s, c) in enumerate(zip(gsamp, gct))}

samples = meta.obs_names.tolist()
cell_types = sorted(set(gct.tolist()))

# --- cell-type proportions obsm (from the meta obs cell_type_* columns) ---
prop_cols = [c for c in meta.obs.columns if c.startswith("cell_type_")]
meta.obsm["cell_type_proportions"] = meta.obs[prop_cols].astype(float).values
meta.uns["cell_type_proportions_columns"] = [c.replace("cell_type_", "") for c in prop_cols]

# --- per-cell-type pseudobulk obsm ---
ctmap = {}
for ctname in cell_types:
    M = np.full((len(samples), len(genes)), np.nan, dtype=np.float32)
    for si, s in enumerate(samples):
        gi = group_of.get((s, ctname))
        if gi is not None:
            M[si] = PB[gi]
    key = "pseudobulk_" + re.sub(r"[^0-9a-zA-Z]+", "_", ctname).strip("_")
    meta.obsm[key] = M
    ctmap[key] = ctname

meta.uns["pseudobulk_genes"] = genes
meta.uns["pseudobulk_celltype_map"] = ctmap
meta.uns["pseudobulk_note"] = (
    "obsm['pseudobulk_<ct>']: per-sample mean of X_log_norm_counts over that "
    "sample's cells of cell type <ct> (HVG genes = uns['pseudobulk_genes']); "
    "NaN if the sample has 0 such cells. obsm['cell_type_proportions'] columns = "
    "uns['cell_type_proportions_columns'] (percent).")

out = os.path.join(REPO, "data", ds, f"{ds}_meta_with_pseudobulk.h5ad")
meta.write(out, compression="gzip")
print(f"{ds}: {meta.shape} | {len(cell_types)} cell types | {len(genes)} genes")
print("  added obsm: cell_type_proportions +", len(ctmap), "pseudobulk_<ct>")
print("  saved", out, f"({os.path.getsize(out)/1e6:.1f} MB)")
