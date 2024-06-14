import random

import numpy as np
import pandas as pd
import patient_representation as pr
import scanpy as sc 
import os
import time
import traceback

start_time = time.time()


random.seed(42)
np.random.seed(42)

## VIASH START
par = {
    "input": "combat_processed.h5ad",
    "output": "combat_represent.h5ad",
    "cell_type_key":"Annotation_major_subset",
    "sample_key":"scRNASeq_sample_ID",
    "batch_effect_covariate":"Pool_ID",
    "batch_covariates":["scRNASeq_sample_ID","Pool_ID"],
    "output_compression": "gzip",
}
## VIASH END

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
        if isinstance(method_instance, pr.tl.CellTypePseudobulk):
            # Set stricter thresholds for sample and cluster size
            method_instance.prepare_anndata(adata, sample_size_threshold=500, cluster_size_threshold=0)
        else:
            method_instance.prepare_anndata(adata, sample_size_threshold=0, cluster_size_threshold=0)
        
        print("Calculating distances")
        method_instance.calculate_distance_matrix(force=True)
        print("Saving results")
        adata = save_results(adata, method_instance, method_name)
        adata.write(RESULT_PATH)

        print("Success", end="\n\n")

    except Exception as e:
        print(f"An error occurred during {method_name}:")
        traceback.print_exc() 
        print("Failed", end="\n\n")

    return adata


RESULT_PATH = par["output"]
ADATA_PATH = par["input"]
CELL_TYPE_KEY = par["cell_type_key"]
SAMPLE_KEY = par["sample_key"]
BATCH_COVARIATES = par["batch_covariates"]
BATCH_EFFECT = par["batch_effect_covariate"]
OUTPUT_NAME = os.path.splitext(os.path.basename(RESULT_PATH))[0]


adata = sc.read(ADATA_PATH)
print("AT START: ")
print(adata)

# Run scPoli manually to save its cell representation
try:
    print("Setting up scPoli")
    scpoli = pr.tl.SCPoli(sample_key=SAMPLE_KEY, cells_type_key=CELL_TYPE_KEY, layer="X_raw_counts")

    print("Preparing adata")
    scpoli.prepare_anndata(adata, sample_size_threshold=0, cluster_size_threshold=0, optimize_adata=False)

    print("Calculating distances")
    scpoli.calculate_distance_matrix(force=True)

    print("Saving results")
    adata = save_results(adata, scpoli, "scpoli")

    print("Saving scPoli cell representation")
    adata.obsm["X_scpoli"] = scpoli.model.get_latent(
        scpoli.adata,
        mean=True
    )
    adata.write(RESULT_PATH)
    print("AFTER SCPOLI: ")
    print(adata)
    print("Success", sep="\n\n")
except Exception as e:
    print(e)
    print("Failed", sep="\n\n")


#TODO:PILOT fails on raw count
#TODO: update and fixs MRVI and add it to methods 

base_layers = [
    ("X_raw_counts", "raw_counts"),
    ("X_pca", "pca"),
    ("X_harmony", "harmony"),
    ("X_scpoli", "scpoli"),
]
methods = [
    (pr.tl.RandomVector, "random_vec", {}),
    (pr.tl.CellTypesComposition, "composition", {})
]
for layer, method_name_suffix in base_layers:
    methods.extend([
        (pr.tl.CellTypePseudobulk, f"ct_pseudobulk_{method_name_suffix}", {"layer": layer}),
        (pr.tl.TotalPseudobulk, f"pseudobulk_{method_name_suffix}", {"layer": layer}),
        (pr.tl.PILOT, f"pilot_{method_name_suffix}", {"patient_state_col": "Outcome", "layer": layer})
    ])
    if layer == "X_scpoli":
        methods.append((pr.tl.WassersteinTSNE, f"wasserstein_{method_name_suffix}", {"replicate_key": CELL_TYPE_KEY, "layer": layer}))

scvi_tools_models = ['scVI', 'scANVI']
for covariate in BATCH_COVARIATES:
    for model in scvi_tools_models:
        layer_name = f"X_{model}_{covariate}"
        methods.extend([
            (pr.tl.TotalPseudobulk, f"pseudobulk_{model}_{covariate}", {"layer": layer_name}),
            (pr.tl.CellTypePseudobulk, f"ct_pseudobulk_{model}_{covariate}", {"layer": layer_name}),
            (pr.tl.WassersteinTSNE, f"wasserstein_{model}_{covariate}", {"replicate_key": CELL_TYPE_KEY, "layer": layer_name}),
            (pr.tl.PILOT, f"pilot_{model}_{covariate}", {"patient_state_col": "Outcome", "layer": layer_name})
        ])

methods.sort(key=lambda x: (x[0].__name__, x[1]))
print("METHODS: ", methods)

for method_class, method_name, kwargs in methods:
    adata = get_representation(adata, method_class, method_name, sample_key=SAMPLE_KEY, cells_type_key=CELL_TYPE_KEY, **kwargs)

print(adata)
# print("Saving layers for GloScope")

# for layer, feature_name in [("X_pca", "PC"), ("X_scVI_Pool_ID", "scvi"), ("X_scANVI_Pool_ID", "scanvi"), ("X_scpoli", "scpoli")]:
#     print("Working with layer", layer)
#     try: 
#         pd.DataFrame(
#             adata.obsm[layer],
#             index=adata.obs_names,
#             columns=[f"{feature_name}{i}" for i in range(adata.obsm[layer].shape[1])]
#         ).to_csv(f"../data/gloscope_input/combat_{layer}.csv")
#     except Exception as e:
#         print("Failed", e)

# print("Saving samples")
# pd.DataFrame(
#     adata.obs[SAMPLE_KEY]
# ).to_csv("../data/gloscope_input/combat_samples.csv", index=False)

#TODO: fix error on phate lib
for cells_per_sample in (200, 500, 700):
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
    adata_subset.write(f"{OUTPUT_NAME}_{cells_per_sample}_cells_per_sample.h5ad")
    print("Trying to run PhEMD")
    adata_subset = get_representation(adata_subset, pr.tl.PhEMD, f"phemd_{cells_per_sample}", sample_key=SAMPLE_KEY, cells_type_key=CELL_TYPE_KEY)
    print("Calculated representation, saving")
    adata_subset.write(f"{OUTPUT_NAME}_{cells_per_sample}_cells_per_sample.h5ad")

print("Done")


end_time = time.time()  # End timing
elapsed_time = end_time - start_time  # Calculate elapsed time
print(f"Elapsed time 10% : {elapsed_time:.2f} seconds")