"""Stratified subsample of a cleaned h5ad for fast smoke tests.

Reads the input in backed mode and pulls only the selected rows into
memory, so it works on files far larger than RAM. Samples up to
--per-sample cells from each group in --sample-key (keeps every sample
represented so the sample/batch keys are still exercised downstream).

Usage:
    python scripts/subsample_h5ad.py \
        --input  data/imm_of_aging/imm_of_aging_cleaned.h5ad \
        --output data/imm_of_aging/imm_of_aging_smoke.h5ad \
        --sample-key sample.sampleKitGuid \
        --per-sample 1300
"""
from __future__ import annotations

import argparse

import anndata as ad
import numpy as np


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--sample-key", required=True)
    p.add_argument("--per-sample", type=int, default=1300)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)

    print(f"Opening {args.input} (backed)")
    a = ad.read_h5ad(args.input, backed="r")
    sample = a.obs[args.sample_key].astype(str).to_numpy()
    print(f"n_obs={a.n_obs:,}  samples={len(set(sample))}")

    keep = []
    for s in sorted(set(sample)):
        idx = np.flatnonzero(sample == s)
        if idx.size > args.per_sample:
            idx = rng.choice(idx, args.per_sample, replace=False)
        keep.append(idx)
    keep = np.sort(np.concatenate(keep))
    print(f"Selected {keep.size:,} cells")

    sub = a[keep].to_memory()
    a.file.close()
    print(f"Subsampled AnnData: {sub.shape}")

    print(f"Writing {args.output}")
    sub.write(args.output)
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
