import random

import numpy as np
import pandas as pd
import patient_representation as pr
import scanpy as sc 

random.seed(42)
np.random.seed(42)


def save_results(adata, method, method_name):
    method.embed(method="UMAP") 
    adata.uns[method_name + "_" + "distances"] = method.adata.uns[method.DISTANCES_UNS_KEY]
    adata.uns[method_name + "_" + "UMAP"] = np.array(method.embeddings["UMAP"])
    adata.uns[method_name + "_" + "samples"] = method.samples 
    return adata

def get_representation(adata, method_class, method_name, **kwargs):
    print(f"Running {method_name}")
    try:
        print("Initializing method")
        method_instance = method_class(**kwargs)

        print("Preparing adata")
        if method_name == "ct_pseudobulk":
            # Set stricter thresholds for sample and cluster size
            method_instance.prepare_anndata(adata, sample_size_threshold=500, cluster_size_threshold=5)
        else:
            method_instance.prepare_anndata(adata, sample_size_threshold=0, cluster_size_threshold=0)
        
        print("Calculating distances")
        method_instance.calculate_distance_matrix(force=True)
        print("Saving results")
        adata = save_results(adata, method_instance, method_name)
        adata.write(RESULT_PATH)

        print("Success", end="\n\n")

    except Exception as e:
        print(e)
        print("Failed", end="\n\n")

    return adata


ADATA_PATH = "../data/onek1k_processed.h5ad"
RESULT_PATH = "../data/onek1k_with_representations.h5ad"
SAMPLE_KEY = "donor_id"
CELL_TYPE_KEY = "cell_type"

adata = sc.read(ADATA_PATH)
print(adata)

# Run scPoli manually to save its cell representation
try:
    print("Setting up scPoli")
    scpoli = pr.tl.SCPoli(sample_key=SAMPLE_KEY, cells_type_key=CELL_TYPE_KEY, layer="raw")

    print("Preparing adata")
    scpoli.prepare_anndata(adata, sample_size_threshold=0, cluster_size_threshold=0)

    print("Calculating distances")
    scpoli.calculate_distance_matrix(force=True)

    print("Saving results")
    adata = save_results(adata, scpoli, "scpoli")

    print("Saving scPoli cell representation")
    adata.obsm["X_scpoli"] = scpoli.model.get_latent(
        adata.X,
        adata.obs[SAMPLE_KEY].values,
        mean=True
    )
    adata.write(RESULT_PATH)
    print("Success", sep="\n\n")
except Exception as e:
    print(e)
    print("Failed", sep="\n\n")

# method class, method name, additional arguments
methods = [
    (pr.tl.RandomVector, "random_vec", {}),
    (pr.tl.TotalPseudobulk, "pseudobulk", {"layer": "raw"}),
    (pr.tl.TotalPseudobulk, "pseudobulk_pca", {"layer": "X_pca"}),
    (pr.tl.TotalPseudobulk, "pseudobulk_harmony", {"layer": "X_harmony"}),
    (pr.tl.CellTypePseudobulk, "ct_pseudobulk", {}),
    (pr.tl.CellTypesComposition, "composition", {}),
    (pr.tl.TotalPseudobulk, "pseudobulk_scvi", {"layer": "X_scVI"}),
    (pr.tl.TotalPseudobulk, "pseudobulk_scanvi", {"layer": "X_scANVI"}),
    (pr.tl.TotalPseudobulk, "pseudobulk_scpoli", {"layer": "X_scpoli"}),
    (pr.tl.WassersteinTSNE, "wasserstein_scvi", {"replicate_key": CELL_TYPE_KEY, "layer": "X_scVI"}),
    (pr.tl.WassersteinTSNE, "wasserstein_scanvi", {"replicate_key": CELL_TYPE_KEY, "layer": "X_scANVI"}),
    (pr.tl.MrVI, "mrvi", {"categorical_nuisance_keys": ["pool_number"], "layer": "raw", "max_epochs": 50}),
    (pr.tl.PILOT, "pilot", {"patient_state_col": "sex", "layer": "X_pca"}),
    (pr.tl.PILOT, "pilot_harmony", {"patient_state_col": "sex", "layer": "X_harmony"}),
    (pr.tl.PILOT, "pilot_scvi", {"patient_state_col": "sex", "layer": "X_scVI"}),
    (pr.tl.PILOT, "pilot_scanvi", {"patient_state_col": "sex", "layer": "X_scANVI"}),
    (pr.tl.PILOT, "pilot_scpoli", {"patient_state_col": "sex", "layer": "X_scpoli"}),
]

for method_class, method_name, kwargs in methods:
    adata = get_representation(adata, method_class, method_name, sample_key=SAMPLE_KEY, cells_type_key=CELL_TYPE_KEY, **kwargs)

print("Saving layers for GloScope")

for layer, feature_name in [("X_pca", "PC"), ("X_harmony", "harmony"), ("X_scVI", "scvi"), ("X_scANVI", "scanvi"), ("X_scpoli", "scpoli")]:
    print("Working with layer", layer)
    try: 
        pd.DataFrame(
            adata.obsm[layer],
            index=adata.obs_names,
            columns=[f"{feature_name}{i}" for i in range(adata.obsm[layer].shape[1])]
        ).to_csv(f"../data/gloscope_input/onek1k_{layer}.csv")
    except Exception as e:
        print("Failed", e)

print("Saving samples")
pd.DataFrame(
    adata.obs[SAMPLE_KEY]
).to_csv("../data/gloscope_input/onek1k_samples.csv", index=False)

for cells_per_sample in (200, 500):
    print("Subsetting cells, trying", cells_per_sample, "cells per sample")
    adata_subset = pr.pp.subsample(
        adata,
        obs_category_col=SAMPLE_KEY,
        min_samples_per_category=cells_per_sample,
        n_obs=cells_per_sample
    #     fraction=0.1
    )
    print("Subset size:", adata_subset.shape)
    print("Saving")
    adata_subset.write(f"../data/onek1k_{cells_per_sample}_cells_per_sample.h5ad")
    print("Trying to run PhEMD")
    adata_subset = get_representation(adata_subset, pr.tl.PhEMD, f"phemd_{cells_per_sample}", sample_key=SAMPLE_KEY, cells_type_key=CELL_TYPE_KEY)
    print("Calculated representation, saving")
    adata_subset.write(f"../data/onek1k_{cells_per_sample}_cells_per_sample.h5ad")

print("Done")
