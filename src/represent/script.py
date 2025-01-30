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
    "input": "data/combat_processed.h5ad",
    "output": "data/combat_represent.h5ad",
    "cell_type_key":"Annotation_major_subset",
    "sample_key":"scRNASeq_sample_ID",
    "batch_covariates":["scRNASeq_sample_ID","Pool_ID","Institute"],
    "celltype_pseudobulk_sample_size_threshold": "0",
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
        if isinstance(method_instance, pr.tl.GroupedPseudobulk):
            adata = pr.pp.filter_small_samples(adata, sample_key=par["sample_key"], sample_size_threshold=int(par["celltype_pseudobulk_sample_size_threshold"]))
        
        method_instance.prepare_anndata(adata)
        
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
OUTPUT_NAME = os.path.splitext(os.path.basename(RESULT_PATH))[0]

adata = sc.read(ADATA_PATH)
print("AT START: ")
print(adata)

# Run scPoli manually to save its cell representation
try:
    print("Setting up scPoli")
    scpoli = pr.tl.SCPoli(sample_key=SAMPLE_KEY, cell_group_key=CELL_TYPE_KEY, layer="X_raw_counts")

    print("Preparing adata")
    scpoli.prepare_anndata(adata, optimize_adata=False)

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

#TODO: update and fixs MRVI and add it to methods 
base_layers = [
    ("X_raw_counts", "raw_counts"),
    ("X_pca", "pca"),
    ("X_harmony", "harmony"),
    ("X_scpoli", "scpoli"),
]
methods = [
    (pr.tl.RandomVector, "random_vec", {}),
    (pr.tl.CellGroupComposition, "composition", {})
]
for layer, method_name_suffix in base_layers:
    methods.extend([
        (pr.tl.GroupedPseudobulk, f"ct_pseudobulk_{method_name_suffix}", {"layer": layer}),
        (pr.tl.Pseudobulk, f"pseudobulk_{method_name_suffix}", {"layer": layer}),
    ])
    if layer != "X_raw_counts":
        methods.append((pr.tl.PILOT, f"pilot_{method_name_suffix}", {"patient_state_col": SAMPLE_KEY, "layer": layer}))
    if layer == "X_scpoli":
        methods.append((pr.tl.WassersteinTSNE, f"wasserstein_{method_name_suffix}", {"replicate_key": CELL_TYPE_KEY, "layer": layer}))

scvi_tools_models = ['scVI', 'scANVI']
for covariate in BATCH_COVARIATES:
    for model in scvi_tools_models:
        layer_name = f"X_{model}_{covariate}"
        methods.extend([
            (pr.tl.Pseudobulk, f"pseudobulk_{model}_{covariate}", {"layer": layer_name}),
            (pr.tl.GroupedPseudobulk, f"ct_pseudobulk_{model}_{covariate}", {"layer": layer_name}),
            (pr.tl.WassersteinTSNE, f"wasserstein_{model}_{covariate}", {"replicate_key": CELL_TYPE_KEY, "layer": layer_name}),
            (pr.tl.PILOT, f"pilot_{model}_{covariate}", {"patient_state_col": SAMPLE_KEY, "layer": layer_name})
        ])

methods.sort(key=lambda x: (x[0].__name__, x[1]))
print("METHODS: ", methods)

for method_class, method_name, kwargs in methods:
    adata = get_representation(adata, method_class, method_name, sample_key=SAMPLE_KEY, cell_group_key=CELL_TYPE_KEY, **kwargs)

print(adata)
print("Done")


end_time = time.time()  # End timing
elapsed_time = end_time - start_time  # Calculate elapsed time
print(f"Elapsed time 10% : {elapsed_time:.2f} seconds")
