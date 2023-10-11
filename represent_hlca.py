import random

import numpy as np
import patient_representation as pr
import pandas as pd
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


ADATA_PATH = "../data/HLCA_processed.h5ad"
RESULT_PATH = "../data/HLCA_with_representations.h5ad"
SAMPLE_KEY = "donor_id"
CELL_TYPE_KEY = "cell_type"

adata = sc.read(ADATA_PATH)
print(adata)

# Run scPoli manually to save its cell representation
try:
    print("Setting up scPoli")
    scpoli = pr.tl.SCPoli(sample_key=SAMPLE_KEY, cells_type_key=CELL_TYPE_KEY, layer="X")

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
    (pr.tl.TotalPseudobulk, "pseudobulk", {"layer": "X_pca"}),
    (pr.tl.TotalPseudobulk, "pseudobulk_pca", {"layer": "X_pca"}),
    (pr.tl.CellTypePseudobulk, "ct_pseudobulk", {}),
    (pr.tl.CellTypesComposition, "composition", {}),
    (pr.tl.TotalPseudobulk, "pseudobulk_scanvi", {"layer": "X_scanvi_emb"}),
    (pr.tl.TotalPseudobulk, "pseudobulk_scpoli", {"layer": "X_scpoli"}),
    (pr.tl.WassersteinTSNE, "wasserstein_scanvi", {"replicate_key": CELL_TYPE_KEY, "layer": "X_scanvi_emb"}),
    (pr.tl.MrVI, "mrvi", {"categorical_nuisance_keys": ["dataset"], "layer": "X", "max_epochs": 100}),
    (pr.tl.PILOT, "pilot", {"patient_state_col": "disease", "layer": "X_pca"})
]

for method_class, method_name, kwargs in methods:
    adata = get_representation(adata, method_class, method_name, sample_key=SAMPLE_KEY, cells_type_key=CELL_TYPE_KEY, **kwargs)

print("Saving layers for GloScope")

for layer, feature_name in [("X_pca", "PC"), ("X_scanvi_emb", "scanvi"), ("X_scpoli", "scpoli")]:
    print("Working with layer", layer)
    try: 
        pd.DataFrame(
            adata.obsm[layer],
            index=adata.obs_names,
            columns=[f"{feature_name}{i}" for i in range(adata.obsm[layer].shape[1])]
        ).to_csv(f"../data/hlca_{layer}.csv")
    except Exception as e:
        print("Failed", e)

print("Saving samples")
pd.DataFrame(
    adata.obs[SAMPLE_KEY]
).to_csv("../data/hlca_samples.csv", index=False)

print("Done")
