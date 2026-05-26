import random

import scanpy as sc
import pandas as pd
import patpy as pr
import numpy as np
from scipy.spatial.distance import pdist, squareform

random.seed(42)
np.random.seed(42)

## VIASH START
par = {
    "input": "data/combat/combat_processed.h5ad",
    "output": "data/combat/combat_sampleclr_pca.csv",
    "sample_key": "scRNASeq_sample_ID",
    "cell_type_key": "Annotation_major_subset",
    "layer": "X_pca",
    "batch_key": "Pool_ID",
    "use_batch_aware_sampler": False,
    "output_dim": 128,
    "n_cells_per_sample": 1000,
    "num_epochs_pretrain": 100,
    "batch_size": 32,
    "device": "cuda",
    "seed": 42,
}
## VIASH END

print("Reading adata")
adata = sc.read(par["input"])
print(adata)

if par["layer"] not in adata.obsm:
    raise ValueError(
        f"Layer '{par['layer']}' not found in adata.obsm "
        f"(available: {list(adata.obsm.keys())})"
    )

batch_key = par["batch_key"] if par["use_batch_aware_sampler"] else None
if par["use_batch_aware_sampler"] and not batch_key:
    raise ValueError("--use_batch_aware_sampler=true requires --batch_key to be set.")

print(
    f"Setting up SampleCLR (self-supervised, layer={par['layer']}, "
    f"batch_aware={par['use_batch_aware_sampler']}, batch_key={batch_key})"
)
model = pr.tl.SampleCLR(
    sample_key=par["sample_key"],
    label_keys=[],
    tasks=[],
    layer=par["layer"],
    output_dim=par["output_dim"],
    n_cells_per_sample=par["n_cells_per_sample"],
    num_epochs_pretrain=par["num_epochs_pretrain"],
    batch_size=par["batch_size"],
    batch_key=batch_key,
    device=par["device"],
    seed=par["seed"],
)

print(f"Training (pretrain only, {par['num_epochs_pretrain']} epochs)")
model.prepare_anndata(
    adata,
    pretrain=True,
    fine_tune=False,
    verbose=True,
)

print("Extracting sample representations")
embeddings = model.get_sample_representations()
print(f"Embeddings shape: {embeddings.shape}")

# Match the contract every other representation component obeys: emit a
# sample × sample distance matrix indexed by sample id so
# aggregate_representations can call ep.pp.neighbors with metric=precomputed
# uniformly. Cosine distance is the natural fit for a contrastive-learning
# embedding — the SampleCLR loss optimises angular separation between
# samples, so |1 - cosine_similarity| matches the geometry the embedding was
# trained for.
print("Computing pairwise cosine distances between sample representations")
distances = squareform(pdist(embeddings.to_numpy(), metric="cosine"))
distance_df = pd.DataFrame(
    distances, index=embeddings.index, columns=embeddings.index
)

print(f"Saving distance matrix ({distance_df.shape}) to {par['output']}")
distance_df.to_csv(par["output"], index=True)
