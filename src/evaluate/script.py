import random
import json
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ehrapy as ep
from scipy import stats
import patient_representation as pr


## VIASH START
par = {
    "input": "data/combat/combat_representations.h5ad",
    "output_dir": "data/combat/results",
    "benchmark_schema": "data/combat/benchmark_schema.json",
    "root_sample": "H00052-Ha001E-PBGa",  # the youngest healthy non-smoker
    "trajectory_variable": "Outcome",
    "inverse_trajectory": True,
    "figure_format": "png",
}
## VIASH END

random.seed(42)
np.random.seed(42)

def get_col_from_adata(adata, col):
    if col in adata.obs.columns:
        return adata.obs[col]
    else:
        return pd.Series(adata[:, col].X.toarray().flatten(), index=adata.obs_names)

meta_adata = ep.io.read_h5ad(par["input"])

output_dir = Path(par["output_dir"])
output_dir.mkdir(parents=True, exist_ok=True)

meta_adata.uns["iroot"] = np.flatnonzero(meta_adata.obs_names == par["root_sample"])[0]

representations = meta_adata.uns["sample_representations"]
print("Number of representations: ", len(representations))
print(representations)

if par["trajectory_variable"] is not None:
    for representation in representations:
        try:
            print(f"Computing diffmap for {representation}")
            ep.tl.diffmap(meta_adata, neighbors_key=f"{representation}_neighbors")
            meta_adata.obsm[f"X_{representation}_diffmap"] = meta_adata.obsm["X_diffmap"]
            ep.tl.dpt(meta_adata, neighbors_key=f"{representation}_neighbors")
            meta_adata.obs.rename(columns={"dpt_pseudotime": f"{representation}_dpt_pseudotime"}, inplace=True)
        except Exception as e:
            print(f"Error computing diffmap for {representation}: {e}")
            meta_adata.obs[f"{representation}_dpt_pseudotime"] = np.zeros(len(meta_adata.obs))
            continue

    trajectory_correlations = []

    for representation in representations:
        target = get_col_from_adata(meta_adata, par["trajectory_variable"])
        
        corr, _ = stats.spearmanr(target, meta_adata.obs[f"{representation}_dpt_pseudotime"], nan_policy="omit")

        if par["inverse_trajectory"]:
            corr = -corr
        trajectory_correlations.append(corr)

    trajectory_metric_df = pd.DataFrame(trajectory_correlations, index=representations, columns=["correlation"])
    trajectory_metric_df.sort_values("correlation", ascending=False, inplace=True)
    trajectory_metric_df.to_csv(output_dir / "trajectory_metric.csv")

    plt.figure(figsize=(5, 10))
    sns.barplot(x="correlation", y=trajectory_metric_df.index, data=trajectory_metric_df)
    plt.tight_layout()
    plt.savefig(output_dir / "trajectory_metric.png", format=par["figure_format"])
    plt.close()

benchmark_schema = json.load(open(par["benchmark_schema"]))
print("Benchmark schema:")
print(benchmark_schema)

results = []

for representation in representations: 
    for covariate_type in benchmark_schema:
        for col in benchmark_schema[covariate_type]:
            task = benchmark_schema[covariate_type][col]
            try:
                if representation == "ehrapy":
                    distances = meta_adata.obsp[f"{representation}_neighbors_distances"].toarray()
                    # TODO: check that not calculated distances are not 0s
                else:
                    distances = meta_adata.obsm[f"{representation}_distances"].values

                result = pr.tl.evaluate_representation(
                    distances=distances,
                    target=get_col_from_adata(meta_adata, col),
                    method="knn",
                    task=task,
                    n_neighbors=3
                )
            except Exception as e:
                print("Representation:", representation)
                print("Covariate:", col)
                print("Task:", task)
                print("Error:", e)
                print()
                continue

            result["representation"] = representation
            result["covariate"] = col
            result["covariate_type"] = covariate_type

            # Inverse technical score to interpret them as batch effect removal
            if covariate_type == "technical":
                result["score"] = 1 - result["score"]
            
            if result["metric"] == "spearman_r":
                result["score"] = abs(result["score"])
            
            results.append(result)
            
knn_results = pd.DataFrame(results)

avg_knn_score_df = pd.DataFrame(
    knn_results.groupby(["representation", "covariate_type"])["score"].aggregate(np.mean)
).reset_index()
avg_knn_score_df = avg_knn_score_df.sort_values("score", ascending=False)
avg_knn_score_df.to_csv(output_dir / "avg_knn_score_df.csv")

relevant_features_weight = 2/3

averaged_scores = avg_knn_score_df.pivot(index="representation", columns="covariate_type", values="score").sort_values("relevant", ascending=False)

averaged_scores["total"] = (averaged_scores["relevant"] * relevant_features_weight + (1 - relevant_features_weight) * averaged_scores["technical"])
averaged_scores.sort_values("total", ascending=False, inplace=True)
averaged_scores.to_csv(output_dir / "averaged_scores.csv")

sns.scatterplot(data=averaged_scores, x="relevant", y="technical", hue="contextual")
plt.legend(loc=(1.05, 0))
plt.tight_layout()
plt.savefig(output_dir / "tradeoff.png", format=par["figure_format"])
plt.close()

knn_results_wide = knn_results.pivot(index="representation", columns="covariate", values="score")
knn_results_wide.to_csv(output_dir / "knn_results.csv")

# Order by median score
covariate_order = knn_results.groupby("covariate")["score"].median().sort_values(ascending=False).index

plt.figure(figsize=(10, 16), dpi=100)
sns.boxplot(
    knn_results[knn_results["covariate"].isin(covariate_order)],
    x="score", y="covariate", orient="h",
    order=covariate_order,
    color="lightblue"
)

plt.gca().tick_params(axis='y', labelsize=30)
plt.gca().tick_params(axis='x', labelsize=30)
sns.swarmplot(knn_results[knn_results["covariate"].isin(covariate_order)], x="score", y="covariate", orient="h", order=covariate_order)
plt.tight_layout()
plt.savefig(output_dir / "knn_results_boxplot.png", format=par["figure_format"])
plt.close()
