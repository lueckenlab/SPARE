import random
import os
import re

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
    "n_cells_per_sample_train_low": None,
    "n_cells_per_sample_train_high": None,
    "n_cells_per_sample_inference": None,
    "num_epochs_pretrain": 100,
    "batch_size": 16,
    "lambda_": 5.0,
    "contrastive_loss_temperature": 0.1,
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


def _resolve_n_cells(value, per_donor_counts, label):
    """Accept an int, a string of digits, or a 'pNN' percentile keyword."""
    if value is None:
        return None
    if isinstance(value, (int, np.integer)):
        return int(value)
    s = str(value).strip()
    if re.fullmatch(r"\d+", s):
        return int(s)
    m = re.fullmatch(r"[pP](\d{1,3})", s)
    if not m:
        raise ValueError(
            f"--{label} must be an int or a 'pNN' percentile keyword (got {value!r})"
        )
    q = int(m.group(1))
    if not 0 <= q <= 100:
        raise ValueError(f"--{label} percentile must be in [0, 100], got p{q}")
    return int(np.percentile(per_donor_counts, q))


per_donor_counts = adata.obs[par["sample_key"]].value_counts().to_numpy()
train_low = _resolve_n_cells(
    par.get("n_cells_per_sample_train_low"),
    per_donor_counts,
    "n_cells_per_sample_train_low",
)
train_high = _resolve_n_cells(
    par.get("n_cells_per_sample_train_high"),
    per_donor_counts,
    "n_cells_per_sample_train_high",
)
inference_n = _resolve_n_cells(
    par.get("n_cells_per_sample_inference"),
    per_donor_counts,
    "n_cells_per_sample_inference",
)

if (train_low is None) != (train_high is None):
    raise ValueError(
        "--n_cells_per_sample_train_low and --n_cells_per_sample_train_high "
        "must both be set or both omitted."
    )

if train_low is not None:
    if train_low > train_high:
        raise ValueError(
            f"n_cells_per_sample_train_low ({train_low}) > "
            f"n_cells_per_sample_train_high ({train_high})"
        )
    n_cells_arg = (train_low, train_high)
    print(
        f"Training cells-per-donor range: {n_cells_arg} "
        f"(from per-donor counts min={per_donor_counts.min()} "
        f"max={per_donor_counts.max()})"
    )
else:
    n_cells_arg = par["n_cells_per_sample"]
    print(f"Training cells-per-donor fixed: {n_cells_arg}")


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
    n_cells_per_sample=n_cells_arg,
    num_epochs_pretrain=par["num_epochs_pretrain"],
    batch_size=par["batch_size"],
    batch_key=batch_key,
    device=par["device"],
    seed=par["seed"],
    # Forwarded to sampleclr.models.ContrastiveModel via **model_kwargs.
    # aggregator_activation defaults to 'relu' upstream (don't override).
    # val_ids defaults to None (no validation patients) — don't pass.
    # Auxiliary attention-regularization losses ("diversity, etc.") all set
    # to 0 so the only loss is the main contrastive loss. Upstream defaults
    # are 1e-4 / 0.0 / 0.0 — only attention_sparsity_lambda is actually
    # nonzero by default; the other two are pinned explicitly here so a
    # future upstream default change can't silently re-introduce them.
    contrastive_loss_temperature=par["contrastive_loss_temperature"],
    lambda_=par["lambda_"],
    attention_sparsity_lambda=0.0,
    attention_orthogonality_lambda=0.0,
    attention_entropy_lambda=0.0,
)

print(f"Training (pretrain only, {par['num_epochs_pretrain']} epochs)")
model.prepare_anndata(
    adata,
    pretrain=True,
    fine_tune=False,
    verbose=True,
)

# prepare_anndata internally calls _compute_sample_representations with the
# defaults (cell_selection='random', no subset_size). Re-run extraction when
# the experiment overrides the inference subset size.
extract_kwargs = {}
if inference_n is not None:
    extract_kwargs["subset_size"] = inference_n
if extract_kwargs:
    print(f"Re-extracting donor embeddings with {extract_kwargs}")
    model._compute_sample_representations(**extract_kwargs)

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
