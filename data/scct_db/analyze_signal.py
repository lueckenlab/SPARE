"""Quick sample-representation signal check on the GSE246613 cohort (scCT-DB,
largest study: TNBC + pembrolizumab, 94 samples, 4 scCT-DB datasets).

Two simple representations:
  1. cell-type composition  (sample x cell-type fraction)
  2. pseudobulk             (sample x gene, CPM-log, HVG, PCA)

For each we test whether the representation separates samples by the
sample-level labels that *vary within this cohort*: timepoint (baseline vs
post), response group, patient, and source scCT-DB dataset (batch sanity).
Metric = leave-one-out kNN label purity vs. the random expectation
(sum of squared class frequencies). A real signal => purity >> random.
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

HERE = os.path.dirname(__file__)
DS = ["scCT_dataset_154-I", "scCT_dataset_155-I", "scCT_dataset_156-I", "scCT_dataset_157-I"]

# ---------- metadata ----------
M = pd.read_csv(os.path.join(HERE, "master_sample_info.tsv"), sep="\t", dtype=str)
M = M.rename(columns={"scCT-DB dataset ID": "dataset", "Sample ID": "sample"})
M = M[M["dataset"].isin(DS)].copy()
M["key"] = M["dataset"] + "::" + M["sample"]
def tp(x):
    x = str(x).lower()
    return "baseline" if "baseline" in x or x.startswith("pre") else "post"
M["timepoint"] = M["Sampling time"].map(tp)
meta = M.set_index("key")[["dataset", "Patient ID", "Response group category",
                            "timepoint", "Cancer subtype"]]
meta.columns = ["dataset", "patient", "response", "timepoint", "cancer"]

# ---------- kNN purity ----------
def knn_purity(X, labels, k=10, metric="euclidean"):
    labels = pd.Series(labels).values
    mask = pd.notna(labels)
    nn = NearestNeighbors(n_neighbors=min(k + 1, len(X)), metric=metric).fit(X)
    _, nbr = nn.kneighbors(X)
    agree = []
    for i in range(len(X)):
        if not mask[i]:
            continue
        nb = [j for j in nbr[i, 1:] if mask[j]]
        if nb:
            agree.append(np.mean([labels[j] == labels[i] for j in nb]))
    return float(np.mean(agree))

def random_purity(labels):
    v = pd.Series(labels).dropna().value_counts(normalize=True)
    return float((v ** 2).sum())

def report(name, X, idx, metric="euclidean"):
    md = meta.reindex(idx)
    print(f"\n=== {name}  (n={len(idx)}, dims={X.shape[1]}) ===")
    print(f"{'label':<14}{'kNN purity':>11}{'random':>9}{'lift':>7}{'classes':>8}")
    out = {}
    for col in ["dataset", "patient", "response", "timepoint"]:
        lab = md[col].values
        if pd.Series(lab).nunique() < 2:
            print(f"{col:<14}{'(constant)':>11}")
            continue
        p = knn_purity(X, lab, metric=metric)
        r = random_purity(lab)
        print(f"{col:<14}{p:>11.3f}{r:>9.3f}{p-r:>7.3f}{pd.Series(lab).nunique():>8}")
        out[col] = (p, r)
    return out

# ---------- 1. composition ----------
C = pd.read_csv(os.path.join(HERE, "composition_with_meta.csv.gz"))
C["key"] = C["dataset"] + "::" + C["sample"]
C = C[C["dataset"].isin(DS)].set_index("key")
ctcols = [c for c in C.columns if c not in
          ("dataset", "sample", "n_cells", "Cancer subtype", "Response group category",
           "Drug type", "Sampling location category", "Sampling time", "Source dataset",
           "Technology", "timepoint")]
Ccomp = C[ctcols]
Ccomp = Ccomp.loc[:, Ccomp.sum() > 0]
comp_res = report("Composition (Bray-Curtis kNN)", Ccomp.values, Ccomp.index, metric="braycurtis")

# ---------- 2. pseudobulk ----------
pb_path = os.path.join(HERE, "gse246613_pseudobulk.csv.gz")
pca_emb = None
if os.path.exists(pb_path):
    PB = pd.read_csv(pb_path, index_col=0)
    # CPM + log1p
    libsize = PB.sum(1).replace(0, np.nan)
    cpm = (PB.T / libsize).T * 1e4
    logcpm = np.log1p(cpm.fillna(0))
    # HVG: top 2000 by variance
    hvg = logcpm.var(0).sort_values(ascending=False).head(2000).index
    Z = StandardScaler().fit_transform(logcpm[hvg].values)
    npc = min(30, Z.shape[0] - 1)
    pca_emb = PCA(n_components=npc, random_state=42).fit_transform(Z)
    pb_res = report("Pseudobulk (CPM-log, 2000 HVG, 30 PC)", pca_emb, list(logcpm.index))
else:
    print("\n[pseudobulk not built yet -- run build_pseudobulk.py]")
    pb_res = {}

# ---------- figure ----------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
# panel A: composition PCA colored by dataset/timepoint
comp_pca = PCA(n_components=2, random_state=42).fit_transform(
    StandardScaler().fit_transform(Ccomp.values))
md = meta.reindex(Ccomp.index)
for grp, sub in pd.DataFrame({"x": comp_pca[:, 0], "y": comp_pca[:, 1],
                              "tp": md["timepoint"].values}).groupby("tp"):
    axes[0].scatter(sub.x, sub.y, s=18, alpha=.8, label=grp)
axes[0].set_title("Composition PCA — by timepoint")
axes[0].set_xlabel("PC1"); axes[0].set_ylabel("PC2"); axes[0].legend(frameon=False, fontsize=8)

# panel B: pseudobulk PCA colored by response (if available)
if pca_emb is not None:
    md2 = meta.reindex(list(logcpm.index))
    dfp = pd.DataFrame({"x": pca_emb[:, 0], "y": pca_emb[:, 1],
                        "resp": md2["response"].fillna("NA").values,
                        "tp": md2["timepoint"].values})
    for grp, sub in dfp.groupby("resp"):
        axes[1].scatter(sub.x, sub.y, s=18, alpha=.8, label=grp)
    axes[1].set_title("Pseudobulk PCA — by response")
    axes[1].set_xlabel("PC1"); axes[1].set_ylabel("PC2"); axes[1].legend(frameon=False, fontsize=8)

# panel C: signal summary bars
labels = ["dataset", "patient", "response", "timepoint"]
xc = np.arange(len(labels)); w = 0.35
comp_lift = [comp_res.get(l, (np.nan, np.nan)) for l in labels]
axes[2].bar(xc - w/2, [c[0] for c in comp_lift], w, label="composition", color="#4C72B0")
axes[2].bar(xc - w/2, [c[1] for c in comp_lift], w, fill=False, edgecolor="k", linewidth=.8)
if pb_res:
    pb_lift = [pb_res.get(l, (np.nan, np.nan)) for l in labels]
    axes[2].bar(xc + w/2, [c[0] for c in pb_lift], w, label="pseudobulk", color="#C44E52")
    axes[2].bar(xc + w/2, [c[1] for c in pb_lift], w, fill=False, edgecolor="k", linewidth=.8)
axes[2].set_xticks(xc); axes[2].set_xticklabels(labels)
axes[2].set_ylim(0, 1); axes[2].set_ylabel("kNN purity (bars=obs, outline=random)")
axes[2].set_title("Signal by label"); axes[2].legend(frameon=False, fontsize=8)
for ax in axes:
    sns.despine(ax=ax); ax.grid(False)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "gse246613_signal.png"), dpi=140)
print("\nSAVED gse246613_signal.png")
