"""Biological interpretation of the SampleCLR sample representation (CD45+).

1. Leiden-cluster the samples using the SampleCLR sample x sample distances.
2. Find the cell-type proportions that differ most across clusters
   (Kruskal-Wallis across leiden groups).
3. For the top discriminating cell type, build a sample x gene pseudobulk
   (mean log-norm expression over that cell type's cells per sample) and rank
   genes between the two largest clusters.
4. Relate clusters to response / timepoint and write a summary.
"""
import os, sys, numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, scanpy as sc
from scipy.stats import kruskal, mannwhitneyu
from sklearn.neighbors import kneighbors_graph
import igraph as ig, leidenalg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

REPO = "/ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark"
DS = "scct_db_cd45pos"
VARIANT = "sampleclr_scanvi_sample_id_t05"   # best response signal among sampleclr
RES = 0.6
HERE = os.path.join(REPO, "data", DS)
OUT = os.path.join(REPO, "data", "scct_db")

# ---------- 1. leiden on sample x sample distances ----------
a = ad.read_h5ad(os.path.join(REPO, "output_scct_db", f"{DS}.sample_representation.output.h5ad"))
D = np.asarray(a.obsm[f"{VARIANT}_distances"], dtype=float)
n = D.shape[0]
k = min(10, n - 1)
A = kneighbors_graph(D, n_neighbors=k, metric="precomputed", mode="connectivity", include_self=False)
A = ((A + A.T) > 0)
src, tgt = A.nonzero()
g = ig.Graph(n=n, edges=list({tuple(sorted((int(s), int(t)))) for s, t in zip(src, tgt)}), directed=False)
part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition,
                                resolution_parameter=RES, seed=42)
a.obs["leiden"] = pd.Categorical([str(x) for x in part.membership])
print("leiden clusters:", a.obs["leiden"].value_counts().to_dict())

# cluster vs metadata
print("\n=== leiden x response ===\n", pd.crosstab(a.obs.leiden, a.obs.response))
print("\n=== leiden x timepoint ===\n", pd.crosstab(a.obs.leiden, a.obs.timepoint))
print("\n=== leiden x dataset ===\n", pd.crosstab(a.obs.leiden, a.obs.dataset))

# ---------- 2. most unequally distributed cell-type proportions ----------
comp_cols = [c for c in a.obs.columns if c.startswith("cell_type_")]
comp = a.obs[comp_cols].astype(float)
comp.columns = [c.replace("cell_type_", "") for c in comp_cols]
groups = a.obs["leiden"].values
rows = []
for ct in comp.columns:
    vals = [comp[ct].values[groups == cl] for cl in sorted(set(groups))]
    vals = [v for v in vals if len(v) > 0]
    try:
        H, p = kruskal(*vals)
    except Exception:
        H, p = np.nan, np.nan
    means = {cl: float(comp[ct].values[groups == cl].mean()) for cl in sorted(set(groups))}
    rows.append({"cell_type": ct, "kruskal_H": H, "p": p, "spread": max(means.values()) - min(means.values()),
                 **{f"mean_cl{cl}": means[cl] for cl in means}})
ct_df = pd.DataFrame(rows).sort_values("kruskal_H", ascending=False)
ct_df.to_csv(os.path.join(OUT, "sampleclr_celltype_by_cluster.csv"), index=False)
print("\n=== Top cell types most unequally distributed across leiden clusters ===")
print(ct_df.head(8).round(3).to_string(index=False))
TOP_CT = ct_df.iloc[0]["cell_type"]

# ---------- 3. cell-type-specific gene expression (sample x celltype pseudobulk) ----------
proc = ad.read_h5ad(os.path.join(HERE, f"{DS}_processed.h5ad"))
proc.X = proc.layers["X_log_norm_counts"]
# subset cells of the top discriminating cell type, pseudobulk per sample (mean log-norm)
sub = proc[proc.obs["cell_type"] == TOP_CT]
print(f"\nTop cell type '{TOP_CT}': {sub.n_obs} cells across {sub.obs['sample_id'].nunique()} samples")
sid = sub.obs["sample_id"].values
samples = sorted(set(sid) & set(a.obs_names))
X = sub.X.tocsr() if sp.issparse(sub.X) else sub.X
pb = np.zeros((len(samples), sub.n_vars))
for i, s in enumerate(samples):
    m = sid == s
    pb[i] = np.asarray(X[m].mean(axis=0)).ravel()
pb = pd.DataFrame(pb, index=samples, columns=sub.var_names)
clmap = a.obs["leiden"].to_dict()
pb_cl = np.array([clmap[s] for s in samples])
# two largest clusters
big = pd.Series(pb_cl).value_counts().index[:2].tolist()
g0 = pb.values[pb_cl == big[0]]; g1 = pb.values[pb_cl == big[1]]
diff = g1.mean(0) - g0.mean(0)
pvals = np.array([mannwhitneyu(g0[:, j], g1[:, j]).pvalue if (g0[:, j].std() + g1[:, j].std()) > 0 else 1.0
                  for j in range(pb.shape[1])])
gene_df = pd.DataFrame({"gene": pb.columns, "logFC_cl{}_minus_cl{}".format(big[1], big[0]): diff, "p": pvals})
gene_df = gene_df.sort_values(by=gene_df.columns[1], key=lambda s: s.abs(), ascending=False)
gene_df.to_csv(os.path.join(OUT, "sampleclr_topct_gene_de.csv"), index=False)
print(f"\n=== {TOP_CT}: genes most different between cluster {big[1]} and {big[0]} ===")
print(gene_df.head(15).round(3).to_string(index=False))

# ---------- 4. figures ----------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
xy = a.obsm[f"{VARIANT}_umap"]
for col, ax, title in [("leiden", axes[0], "leiden cluster"),
                       ("response", axes[1], "response"),
                       ("timepoint", axes[2], "timepoint")]:
    cats = sorted(a.obs[col].astype(str).unique()); cmap = plt.get_cmap("tab10")
    for i, c in enumerate(cats):
        m = a.obs[col].astype(str).values == c
        ax.scatter(xy[m, 0], xy[m, 1], s=40, color=cmap(i % 10), label=c, alpha=.85, linewidths=0)
    ax.set_title(f"SampleCLR sample UMAP — {title}"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "sampleclr_leiden_umap.png"), dpi=140, bbox_inches="tight")

# cell-type proportion boxplots for top 4 cell types
top4 = ct_df.head(4)["cell_type"].tolist()
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
for ax, ct in zip(axes, top4):
    dfp = pd.DataFrame({"leiden": groups, "prop": comp[ct].values})
    sns.boxplot(data=dfp, x="leiden", y="prop", ax=ax, color="#bbbbbb", fliersize=0)
    sns.stripplot(data=dfp, x="leiden", y="prop", ax=ax, color="#333", size=3, alpha=.6)
    ax.set_title(ct[:32], fontsize=9); ax.set_ylim(0, None)
    sns.despine(ax=ax); ax.grid(False)
fig.suptitle("Cell-type proportions by SampleCLR leiden cluster (CD45+)", y=1.03)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "sampleclr_celltype_props.png"), dpi=140, bbox_inches="tight")

# top genes heatmap (top cell type, two clusters)
topg = gene_df.head(20)["gene"].tolist()
hm = pd.DataFrame(pb[topg].values, index=[f"{s[:18]}|cl{clmap[s]}" for s in samples], columns=topg)
hm = hm.loc[[i for i, s in zip(hm.index, samples) if pb_cl[samples.index(s)] in big]]
order = np.argsort([clmap[s] for s in samples if clmap[s] in big])
fig, ax = plt.subplots(figsize=(9, max(5, 0.18 * len(hm))))
sns.heatmap(hm.iloc[order], cmap="viridis", ax=ax, cbar_kws={"label": "log-norm"})
ax.set_title(f"{TOP_CT}: top genes (cl{big[0]} vs cl{big[1]})", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "sampleclr_topct_genes_heatmap.png"), dpi=135, bbox_inches="tight")

print("\nSAVED figures: sampleclr_leiden_umap.png, sampleclr_celltype_props.png, sampleclr_topct_genes_heatmap.png")
print("VARIANT:", VARIANT, "| TOP_CT:", TOP_CT, "| clusters compared:", big)
