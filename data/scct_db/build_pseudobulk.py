"""Stream scCT-DB dense gene x cell Matrix.csv.gz files row-wise and build a
sample x gene pseudobulk (summed counts) for the GSE246613 cohort.

The matrices are dense CSV (genes as rows, cell barcodes as columns), too large
to load densely. We never materialize the full matrix: for each gene row we sum
its per-cell counts into sample buckets using the Cell info barcode->sample map.
"""
import gzip, os, sys, time
import numpy as np
import pandas as pd

RAW = os.path.join(os.path.dirname(__file__), "raw")
DATASETS = ["scCT_dataset_154-I", "scCT_dataset_155-I", "scCT_dataset_156-I", "scCT_dataset_157-I"]

def load_cellmap(ds):
    ci = pd.read_csv(os.path.join(RAW, ds, "Cell_info.txt"), sep="\t", dtype=str)
    # Barcode, Cell type, Sample ID
    return dict(zip(ci["Barcode"], ci["Sample ID"]))

def stream_pseudobulk(ds):
    cellmap = load_cellmap(ds)
    path = os.path.join(RAW, ds, "Matrix.csv.gz")
    t0 = time.time()
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
        cells = header[1:]  # first col is empty (gene name col)
        # sample per matrix-column; -1 if barcode not annotated
        samples = [cellmap.get(c, None) for c in cells]
        uniq = sorted({s for s in samples if s is not None})
        s_index = {s: i for i, s in enumerate(uniq)}
        col_sample = np.array([s_index.get(s, -1) for s in samples], dtype=np.int64)
        valid = col_sample >= 0
        col_valid = col_sample[valid]
        nsamp = len(uniq)
        genes = []
        rows = []  # list of per-gene length-nsamp vectors
        for ln, line in enumerate(fh):
            i = line.find(",")
            gene = line[:i]
            arr = np.fromstring(line[i + 1:], sep=",")
            if arr.shape[0] != len(samples):
                # ragged line guard
                arr = np.resize(arr, len(samples))
            pb = np.bincount(col_valid, weights=arr[valid], minlength=nsamp)
            genes.append(gene)
            rows.append(pb.astype(np.float32))
            if ln % 5000 == 0 and ln:
                print(f"  {ds}: {ln} genes, {time.time()-t0:.0f}s", flush=True)
    mat = np.vstack(rows)  # genes x samples
    df = pd.DataFrame(mat, index=genes, columns=uniq).T  # samples x genes
    df = df.groupby(level=0, axis=1).sum()  # collapse duplicate gene symbols
    print(f"{ds}: {df.shape[0]} samples x {df.shape[1]} genes, {time.time()-t0:.0f}s", flush=True)
    return df

def main():
    parts = []
    for ds in DATASETS:
        df = stream_pseudobulk(ds)
        df.index = [f"{ds}::{s}" for s in df.index]
        parts.append(df)
    # union genes, fill 0
    allg = sorted(set().union(*[set(p.columns) for p in parts]))
    parts = [p.reindex(columns=allg, fill_value=0) for p in parts]
    full = pd.concat(parts, axis=0)
    out = os.path.join(os.path.dirname(__file__), "gse246613_pseudobulk.csv.gz")
    full.to_csv(out, compression="gzip")
    print("SAVED", out, full.shape)

if __name__ == "__main__":
    main()
