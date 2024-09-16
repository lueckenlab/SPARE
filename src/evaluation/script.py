import traceback
import scanpy as sc
from subprocess import run
import numpy as np
import pandas as pd
import sys
import patient_representation as pr
import json
import matplotlib.pyplot as plt
import seaborn as sns
from plottable import ColumnDefinition, Table
from plottable.cmap import normed_cmap
from plottable.plots import bar
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
import ehrapy as ep
import os 

## VIASH START
par = {
    "input": "../../data/combat_represent_200_with_qc_metrics.h5ad",
    "cell_type_key":"Annotation_major_subset",
    "output_compression": "gzip",
    "sample_key":"scRNASeq_sample_ID",
    "metadata_path": "../../data/combat_metadata_200.csv",
    # "benchmark_schema" : {
    # "technical": ["Institute", "Pool_ID"],
    # # "technical": ["Institute", "Pool_ID", "n_cells", "median_QC_ngenes"],
    # "clinical": ["Death28", "Outcome", "Source"],
    # # "biological": ["Annotation_major_subset_CD4"],
    # },
    # "cols_with_tasks" : {
    #     # "Annotation_major_subset_CD4": "regression",
    #     "Institute": "classification",
    #     "Pool_ID": "classification",
    #     # "n_cells": "regression",
    #     # "median_QC_ngenes": "regression",
    #     "Death28": "classification",
    #     "Outcome": "ranking",
    #     "Source": "classification",
    # },
    "benchmark_schema_file": "benchmark_schema.json",
    "cols_with_tasks_file": "cols_with_tasks.json", 
    "samples_metadata_cols": ["Source", "Outcome", "Death28", "Institute", "Pool_ID"],
    "output_base_name": "evaluation_results"
}
## VIASH END

ADATA_PATH = par["input"]
CELL_TYPE_KEY = par["cell_type_key"]
SAMPLE_KEY = par["sample_key"]
SAMPLES_METADATA_COLS = par["samples_metadata_cols"]
OUTPUT_BASE_NAME = par["output_base_name"]
METADATA_PATH = par["metadata_path"]


# Load schema from JSON files
with open(par["benchmark_schema_file"], 'r') as f:
    BENCHMARK_SCHEMA = json.load(f)

with open(par["cols_with_tasks_file"], 'r') as f:
    COLS_WITH_TASKS = json.load(f)

adata = sc.read_h5ad(ADATA_PATH)

print("____________________ADATA FIRST")
print(adata)
print("____________________ADATA FIRST")

# Extract metadata
pat_instance = pr.tl.TotalPseudobulk(sample_key=SAMPLE_KEY, cells_type_key=CELL_TYPE_KEY)
pat_instance.prepare_anndata(adata, sample_size_threshold=0, cluster_size_threshold=0)
metadata = pat_instance._extract_metadata(SAMPLES_METADATA_COLS)

print("____________________metadata ")
print(metadata)
print("____________________metadata ")

print("____________________metadata columns")
print(metadata.columns)
print("____________________metadata columns\n")

figures_dir = "figures"
tables_dir = "tables"
os.makedirs(figures_dir, exist_ok=True)
os.makedirs(tables_dir, exist_ok=True)


# Align Representations
def align_representations(adata, meta_adata, samples, methods, cols_of_interest):
    for method in methods:
        try:
            representation_samples = adata.uns[f"{method}_samples"].tolist()
            samples_order = [representation_samples.index(sample) for sample in samples if sample in representation_samples]
                        
            if not (adata.uns[f"{method}_samples"][samples_order] == samples).all():  
                print(f"Warning: Order of samples is not correct for method {method}.\n")  
                continue 
            
            assert (adata.uns[f"{method}_samples"][samples_order] == samples).all(), "Order of samples is not correct"

            meta_adata.obsm[f"{method}_UMAP"] = adata.uns[f"{method}_UMAP"][samples_order]
            meta_adata.obsm["umap"] = meta_adata.obsm[f"{method}_UMAP"]
            meta_adata.obsm[f"{method}_distances"] = adata.uns[f"{method}_distances"][samples_order][:, samples_order]

            ep.pp.neighbors(meta_adata, use_rep=f"{method}_distances", key_added=f"{method}_neighbors", metric="precomputed")
            ep.tl.leiden(meta_adata, key_added=f"{method}_leiden", neighbors_key=f"{method}_neighbors")

            fig = ep.pl.umap(meta_adata, color=[f"{method}_leiden"] + cols_of_interest, return_fig=True)
            fig.suptitle(method, fontsize=20)
            
            fig.savefig(os.path.join(figures_dir, f"{OUTPUT_BASE_NAME}_{method}_UMAP.png"))
            
        except Exception as e:  
            print(f"An error occurred with method {method}: {e}\n") 
            continue  
    return meta_adata

# Get common samples across methods
representations_methods = ["pseudobulk_pca"]
combat_samples = list(set(adata.uns[f"{method}_samples"]) for method in representations_methods)
combat_samples = list(set.intersection(*combat_samples))
print("Common Combat Samples:\n", combat_samples)

# Align metadata with representations
combat_meta_adata = ep.io.df_to_anndata(metadata.loc[combat_samples])
combat_meta_adata = ep.pp.encode(combat_meta_adata, autodetect=True)
combat_meta_adata = align_representations(adata=adata, meta_adata=combat_meta_adata, samples=combat_samples, methods=representations_methods, cols_of_interest=SAMPLES_METADATA_COLS)

print("__________________________combat_meta_adata FIRST")
print(combat_meta_adata)
print("__________________________combat_meta_adata FIRST\n")

# Evaluate Representations
def evaluate_representations(combat_meta_adata, methods, benchmark_schema, cols_with_tasks):
    results = []

    for method in methods: 
        for covariate_type in benchmark_schema:
            for col in benchmark_schema[covariate_type]:
                task = cols_with_tasks[col]
                print(f"Evaluating Method: {method}, Covariate: {col}, Task: {task}\n") 
                 #later on printing should be done via proper logging               
                try:
                    result = pr.tl.evaluate_representation(
                        distances=combat_meta_adata.obsm[f"{method}_distances"],
                        target=metadata[col],
                        # target=metadata.get(col),
                        method="knn",
                        task=task,
                    )
                    if result is None:
                        print(f"No result for method {method}, covariate {col}.\n") 
                        continue 
                except Exception as e:
                    print(f"Method: {method}, Col: {col}, Task: {task}, Error: {e}\n")
                    traceback.print_exc()
                    continue  

                result["representation"] = method
                result["covariate"] = col
                result["covariate_type"] = covariate_type

                if result["metric"] == "spearman_r":
                    result["score"] = abs(result["score"])

                if covariate_type == "technical":
                    result["score"] = 1 - result["score"]

                results.append(result)

    return pd.DataFrame(results)

# Evaluate all representations
knn_results = evaluate_representations(combat_meta_adata, representations_methods, BENCHMARK_SCHEMA, COLS_WITH_TASKS)
knn_results.sort_values("score", ascending=False, inplace=True)


# Plot Results
def plot_knn_results(knn_results, output_base_name):
    if not knn_results.empty: 
        plt.figure(figsize=(10, 20))
        sns.barplot(data=knn_results, y="covariate", x="score", orient="h", hue="representation", palette="tab20")
        plt.xlim(0, 1.05)
        plt.title("KNN-score", fontsize=24)
        
        plt.savefig(os.path.join(figures_dir, f"{output_base_name}_knn_score.png"))

    else: 
        print("No results to plot.\n")

def plot_results_table(knn_results, benchmark_schema, output_base_name):
    if not knn_results.empty:
        knn_results_wide = knn_results.pivot(index="representation", columns="covariate", values="score")
        cols_order = ["total"]

        for covariate_type in benchmark_schema:
            type_cols = benchmark_schema[covariate_type]
            # Check if all required columns are in the DataFrame before adding them to cols_order
            if all(col in knn_results_wide.columns for col in type_cols):
                knn_results_wide[covariate_type] = knn_results_wide[type_cols].mean(axis=1)
                cols_order.extend(type_cols)
                cols_order.append(covariate_type)
            else:
                missing_cols = [col for col in type_cols if col not in knn_results_wide.columns]
                print(f"Warning: Missing columns {missing_cols} in {covariate_type}\n")

        cmap = LinearSegmentedColormap.from_list(
            name="bugw", colors=["#FF9693", "#f2fbd2", "#c9ecb4", "#93d3ab", "#35b0ab"], N=256
        )

        col_defs = [ColumnDefinition(
            "total", width=0.7, plot_fn=bar,
            plot_kw={"cmap": cmap, "plot_bg_bar": True, "annotate": True, "height": 0.5, "lw": 0.5, "formatter": lambda x: round(x, 2)}
        )]

        for covariate_type in benchmark_schema:
            type_cols = benchmark_schema[covariate_type]
            if all(col in knn_results_wide.columns for col in type_cols):
                for col in type_cols:
                    col_def = ColumnDefinition(
                        name=col, width=0.75, formatter=lambda x: round(x, 2),
                        textprops={"ha": "center", "bbox": {"boxstyle": "circle", "pad": 0.35}},
                        cmap=normed_cmap(knn_results["score"], cmap=plt.cm.PiYG, num_stds=2.5),
                        group=covariate_type
                    )
                    col_defs.append(col_def)

                knn_results_wide[covariate_type] = knn_results_wide[type_cols].mean(axis=1)
                col_defs.append(ColumnDefinition(
                    covariate_type, width=0.7, plot_fn=bar,
                    plot_kw={"cmap": cmap, "plot_bg_bar": True, "annotate": True, "height": 0.5, "lw": 0.5, "formatter": lambda x: round(x, 2)}
                ))

        clin_weight = 2 / 3
        if "clinical" in knn_results_wide.columns and "technical" in knn_results_wide.columns:
            knn_results_wide["total"] = clin_weight * knn_results_wide["clinical"] + (1 - clin_weight) * knn_results_wide["technical"]

        # Ensure cols_order only contains columns that exist in knn_results_wide
        cols_order = [col for col in cols_order if col in knn_results_wide.columns]

        if "total" in knn_results_wide.columns:  # Only proceed if "total" was calculated
            fig, ax = plt.subplots(figsize=(22, 10))
            Table(knn_results_wide[cols_order].sort_values("total", ascending=False),
                  column_definitions=tuple(col_defs), ax=ax)

            plt.savefig(os.path.join(tables_dir, f"{output_base_name}_results_table.png"))

        else:
            print("Cannot plot results table because 'total' could not be calculated.\n")
    else:
        print("No results to plot.\n")



# Save and plot results
plot_knn_results(knn_results, OUTPUT_BASE_NAME)
plot_results_table(knn_results, BENCHMARK_SCHEMA, OUTPUT_BASE_NAME)

# Save knn results to file
knn_results.to_csv(os.path.join(tables_dir, f"{OUTPUT_BASE_NAME}_knn_results.csv"), index=False)