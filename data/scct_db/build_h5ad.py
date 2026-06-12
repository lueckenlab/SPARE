"""Assemble the GSE246613 cohort (4 scCT-DB datasets) into one cleaned
cells x genes AnnData for the patient-representation benchmark.

Each source Matrix.csv.gz is a dense gene x cell counts CSV (genes as rows,
cell barcodes as columns). We stream it row-wise (never densify the whole
matrix), collect non-zeros into a COO, collapse duplicate gene symbols by
summing, and transpose to cells x genes. Per-cell obs (cell type, sample,
patient, timepoint, response, sort) come from `Cell info` + the master sample
table. Output: scct_db_cleaned.h5ad with raw counts in X and layers.
"""
import gzip, os, time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
DS = ["scCT_dataset_154-I", "scCT_dataset_155-I", "scCT_dataset_156-I", "scCT_dataset_157-I"]

# sample-level metadata (response, sort, timepoint, patient) from master table
M = pd.read_csv(os.path.join(HERE, "master_sample_info.tsv"), sep="\t", dtype=str)
M = M[M["scCT-DB dataset ID"].isin(DS)].copy()
M["timepoint"] = M["Sampling time"].str.lower().map(
    lambda x: "baseline" if "baseline" in str(x) else "post")
SAMPLE_META = {}  # (dataset, sample_id) -> dict
for _, r in M.iterrows():
    SAMPLE_META[(r["scCT-DB dataset ID"], r["Sample ID"])] = dict(
        patient=r["Patient ID"], timepoint=r["timepoint"],
        response=r["Response group category"],
        sort=r["Cell sorting"] if pd.notna(r["Cell sorting"]) else "unknown")


def build_one(ds):
    t0 = time.time()
    ci = pd.read_csv(os.path.join(RAW, ds, "Cell_info.txt"), sep="\t", dtype=str)
    bc2sample = dict(zip(ci["Barcode"], ci["Sample ID"]))
    bc2ct = dict(zip(ci["Barcode"], ci["Cell type"]))
    path = os.path.join(RAW, ds, "Matrix.csv.gz")
    with gzip.open(path, "rt") as fh:
        cells = fh.readline().rstrip("\n").split(",")[1:]
        ncell = len(cells)
        gene_idx, genes = {}, []
        rows, cols, data = [], [], []
        for ln, line in enumerate(fh):
            i = line.find(",")
            gene = line[:i]
            gi = gene_idx.get(gene)
            if gi is None:
                gi = len(genes); gene_idx[gene] = gi; genes.append(gene)
            arr = np.fromstring(line[i + 1:], sep=",")
            if arr.shape[0] != ncell:
                arr = np.resize(arr, ncell)
            nz = np.nonzero(arr)[0]
            if nz.size:
                rows.append(np.full(nz.size, gi, dtype=np.int32))
                cols.append(nz.astype(np.int32))
                data.append(arr[nz].astype(np.float32))
            if ln % 8000 == 0 and ln:
                print(f"  {ds}: {ln} gene-rows, {time.time()-t0:.0f}s", flush=True)
    rows = np.concatenate(rows); cols = np.concatenate(cols); data = np.concatenate(data)
    G = sp.coo_matrix((data, (rows, cols)), shape=(len(genes), ncell)).tocsr()
    G.sum_duplicates()
    X = G.T.tocsr()  # cells x genes
    obs = pd.DataFrame(index=[f"{ds}|{b}" for b in cells])
    samp = [bc2sample.get(b, "NA") for b in cells]
    obs["scct_sample"] = samp
    obs["sample_id"] = [f"{ds}|{s}" for s in samp]
    obs["cell_type"] = [bc2ct.get(b, "Not assigned") for b in cells]
    obs["dataset"] = ds
    for fld in ["patient", "timepoint", "response", "sort"]:
        obs[fld] = [SAMPLE_META.get((ds, s), {}).get(fld, "NA") for s in samp]
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    print(f"{ds}: {a.shape} nnz={X.nnz/1e6:.0f}M  {time.time()-t0:.0f}s", flush=True)
    return a


def main():
    parts = [build_one(ds) for ds in DS]
    adata = ad.concat(parts, join="outer", index_unique=None, fill_value=0)
    adata.X = adata.X.tocsr()
    adata.X.eliminate_zeros()
    # standard light QC (mirror ticatlas clean_data)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    adata.obs["patient_id"] = adata.obs["patient"]
    adata.obs["source"] = "GSE246613"
    adata.layers["X_raw_counts"] = adata.X.copy()
    adata.var_names_make_unique()
    print("FINAL", adata)
    print("samples:", adata.obs["sample_id"].nunique(),
          "| response:", adata.obs["response"].value_counts().to_dict(),
          "| sort:", adata.obs["sort"].value_counts().to_dict())
    out = os.path.join(HERE, "scct_db_cleaned.h5ad")
    adata.write(out, compression="gzip")
    print("SAVED", out)


if __name__ == "__main__":
    main()
