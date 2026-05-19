"""Render a self-contained HTML EDA report for an AIFI cohort.

Reads obs directly via h5py (anndata's backed mode is too slow when obs
has tens of millions of rows + many categoricals). Materialises obs as a
pandas DataFrame by decoding category codes manually, then computes:

- shape + sparsity
- candidate sample / cell_type / batch / visit / sex / age / condition columns
- per-donor cell-count distribution
- cell-type composition
- per-donor × visit matrix (if longitudinal)
- a suggested dataset_info.yaml stanza

Plots follow CLAUDE.md style: no grid, despine, normalised metrics 0-1.
Output goes to reports/<id>_eda.html.

Usage:
    python scripts/eda_aifi.py --cohort imm_of_aging
    python scripts/eda_aifi.py --cohort sound_life
"""

from __future__ import annotations

import argparse
import base64
import io
import sys
import time
from html import escape
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parents[1]

COHORTS = {
    "imm_of_aging": REPO_ROOT / "data" / "imm_of_aging" / "imm_of_aging.h5ad",
    "sound_life": REPO_ROOT / "data" / "sound_life" / "sound_life.h5ad",
}

DONOR_CANDIDATES = [
    "subject.subjectGuid", "donor_id", "donor", "subjectGuid",
    "sample_id", "subject_id", "DonorID",
]
SAMPLE_CANDIDATES = [
    "sample.sampleKitGuid", "sampleKitGuid", "pbmc_sample_id",
    "sample_id", "sampleId", "sample",
]
CELLTYPE_CANDIDATES = [
    "AIFI_L3", "AIFI_L2", "AIFI_L1",
    "cell_type", "celltype", "annotation", "cell_type_l3",
]
BATCH_CANDIDATES = [
    "batch_id", "batchID", "batch",
    "pool_id", "chip_id", "well_id",
    "cohortGuid", "cohort", "library_id",
]
VISIT_CANDIDATES = [
    "sample.visitName", "visit", "timepoint", "time_point", "visit_name", "day",
]
SEX_CANDIDATES = ["subject.biologicalSex", "sex", "biological_sex"]
AGE_CANDIDATES = [
    "subject.ageAtFirstDraw", "sample.subjectAgeAtDraw",
    "subject.ageGroup", "age_group", "subject.ageAtEnrollment",
    "age", "ageAtEnrollment",
]
CONDITION_CANDIDATES = [
    "subject.cmv", "cmv_status", "cmv",
    "vaccine", "response", "outcome", "condition",
]


def find_first(cols: list[str], candidates: list[str]) -> str | None:
    for c in candidates:
        if c in cols:
            return c
    return None


def candidate_block(cols: list[str], candidates: list[str]) -> list[str]:
    return [c for c in candidates if c in cols]


def read_obs_h5(path: Path, max_cardinality_for_object: int = 200) -> tuple[pd.DataFrame, dict]:
    """Materialise obs from an .h5ad file via h5py.

    Categorical columns decoded as pandas Categorical (codes + categories).
    Object columns with cardinality > max_cardinality_for_object are skipped
    to avoid spending GBs on per-cell barcode strings.
    """
    meta: dict = {}
    columns: dict[str, pd.Series | pd.Categorical] = {}
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        n = None
        for c in sorted(obs.keys()):
            if c.startswith("_"):
                continue
            node = obs[c]
            if isinstance(node, h5py.Group) and "categories" in node:
                cats = node["categories"][...]
                codes = node["codes"][...]
                if cats.dtype.kind in ("O", "S"):
                    cats = np.array([s.decode() if isinstance(s, bytes) else s for s in cats])
                if cats.size > max_cardinality_for_object and cats.size > 0.5 * codes.size:
                    # Things like per-cell `cell_name` / `original_barcodes`; skip.
                    meta.setdefault("skipped_high_cardinality", []).append((c, int(cats.size)))
                    continue
                columns[c] = pd.Categorical.from_codes(codes, categories=cats)
                if n is None:
                    n = len(codes)
            elif hasattr(node, "shape"):
                arr = node[...]
                if arr.dtype.kind == "O":
                    if arr.size > 0 and isinstance(arr.flat[0], (bytes, bytearray)):
                        arr = np.array([s.decode() if isinstance(s, bytes) else s for s in arr])
                    if arr.size > max_cardinality_for_object * 10:
                        meta.setdefault("skipped_high_cardinality", []).append((c, int(arr.size)))
                        continue
                columns[c] = pd.Series(arr)
                if n is None:
                    n = arr.size
        x = f["X"]
        if "data" in x:
            shape_attr = x.attrs.get("shape", None)
            meta["nnz"] = int(x["data"].size)
            meta["x_dtype"] = str(x["data"].dtype)
        else:
            shape_attr = x.shape
            meta["x_dtype"] = str(x.dtype)
        if shape_attr is not None:
            meta["shape"] = tuple(int(v) for v in shape_attr)
        else:
            meta["shape"] = (n, None)
        meta["uns_keys"] = sorted(f.get("uns", {}).keys()) if "uns" in f else []
        meta["obsm_keys"] = sorted(f.get("obsm", {}).keys()) if "obsm" in f else []
        meta["layers"] = sorted(f.get("layers", {}).keys()) if "layers" in f else []
        meta["var_n"] = int(f["var"]["_index"].size) if "_index" in f["var"] else None
    df = pd.DataFrame(columns)
    return df, meta


def fig_to_b64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def img_tag(fig: plt.Figure, alt: str) -> str:
    return f'<img alt="{escape(alt)}" src="data:image/png;base64,{fig_to_b64(fig)}"/>'


def make_donor_cells_hist(donor_counts: pd.Series) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.hist(donor_counts.values, bins=40, color="steelblue", edgecolor="white")
    ax.set_xlabel("cells per donor")
    ax.set_ylabel("donors")
    ax.set_title(f"Cells per donor (n={donor_counts.size:,})")
    ax.grid(False)
    sns.despine(ax=ax)
    return fig


def make_celltype_bar(ct_counts: pd.Series, top: int = 30) -> plt.Figure:
    ct = ct_counts.head(top)
    fig, ax = plt.subplots(figsize=(7, max(3.5, 0.25 * len(ct))))
    ax.barh(ct.index[::-1], ct.values[::-1], color="cornflowerblue")
    ax.set_xlabel("cells")
    ax.set_title(f"Top {top} cell types (of {ct_counts.size:,})")
    ax.grid(False)
    sns.despine(ax=ax)
    return fig


def make_donor_visit_heatmap(obs: pd.DataFrame, donor_key: str, visit_key: str) -> plt.Figure:
    g = obs.groupby([donor_key, visit_key], observed=True).size().unstack(fill_value=0)
    g = g.loc[g.sum(axis=1).sort_values(ascending=False).head(60).index]
    fig, ax = plt.subplots(figsize=(max(5, 0.35 * g.shape[1]), max(4, 0.12 * g.shape[0])))
    sns.heatmap(np.log10(g + 1), cmap="viridis", ax=ax, cbar_kws={"label": "log10(cells+1)"})
    ax.set_xlabel(visit_key)
    ax.set_ylabel(donor_key)
    ax.set_title("Donor × visit cell count (top 60 donors)")
    return fig


def make_value_count_bar(s: pd.Series, title: str, top: int = 20) -> plt.Figure:
    counts = s.value_counts(dropna=False).head(top)
    fig, ax = plt.subplots(figsize=(6, max(2.5, 0.22 * len(counts))))
    ax.barh(counts.index.astype(str)[::-1], counts.values[::-1], color="slateblue")
    ax.set_xlabel("cells")
    ax.set_title(title)
    ax.grid(False)
    sns.despine(ax=ax)
    return fig


def make_age_hist(s: pd.Series) -> plt.Figure:
    s_num = pd.to_numeric(s, errors="coerce").dropna()
    if s_num.empty:
        return make_value_count_bar(s, "age")
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.hist(s_num.values, bins=30, color="goldenrod", edgecolor="white")
    ax.set_xlabel("age")
    ax.set_ylabel("cells")
    ax.set_title("Age distribution (per cell)")
    ax.grid(False)
    sns.despine(ax=ax)
    return fig


def df_to_html(df: pd.DataFrame, max_rows: int = 60) -> str:
    return df.head(max_rows).to_html(border=0, classes="dt")


def section(title: str, body_html: str) -> str:
    return f"<h2>{escape(title)}</h2>\n{body_html}\n"


def render_html(cohort: str, sections: list[str]) -> str:
    style = """
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; max-width: 1100px; margin: 2em auto; color: #222; }
    h1 { border-bottom: 1px solid #ccc; padding-bottom: .2em; }
    h2 { margin-top: 2em; color: #333; }
    code, pre { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }
    code { background: #f5f5f5; padding: 1px 4px; border-radius: 3px; }
    pre  { background: #f7f7f7; padding: .8em; border-radius: 4px; overflow-x: auto; font-size: 0.85em; }
    table { border-collapse: collapse; margin: .8em 0; font-size: 0.9em; }
    th, td { border: 1px solid #ddd; padding: 4px 8px; text-align: left; vertical-align: top; }
    th { background: #f0f0f0; }
    img { max-width: 100%; height: auto; }
    .candidates { color: #555; }
    """
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{escape(cohort)} EDA</title><style>{style}</style></head>
<body>
<h1>{escape(cohort)} — EDA report</h1>
{''.join(sections)}
</body></html>
"""


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", required=True, choices=list(COHORTS))
    p.add_argument("--max-categorical-values", type=int, default=15)
    args = p.parse_args()

    path = COHORTS[args.cohort]
    if not path.exists():
        print(f"missing {path}", file=sys.stderr)
        return 1
    out_path = REPO_ROOT / "reports" / f"{args.cohort}_eda.html"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"=== Reading obs from {path.name} via h5py ===", flush=True)
    t0 = time.time()
    obs, meta = read_obs_h5(path)
    print(f"obs read: shape={obs.shape} in {time.time()-t0:.1f}s", flush=True)

    n_obs, n_vars = meta["shape"]
    obs_cols = list(obs.columns)
    donor_key = find_first(obs_cols, DONOR_CANDIDATES)
    sample_key = find_first(obs_cols, SAMPLE_CANDIDATES) or donor_key
    celltype_key = find_first(obs_cols, CELLTYPE_CANDIDATES)
    batch_key = find_first(obs_cols, BATCH_CANDIDATES)
    visit_key = find_first(obs_cols, VISIT_CANDIDATES)
    sex_key = find_first(obs_cols, SEX_CANDIDATES)
    age_key = find_first(obs_cols, AGE_CANDIDATES)
    cond_key = find_first(obs_cols, CONDITION_CANDIDATES)

    sections: list[str] = []

    overview = pd.DataFrame({
        "key": ["n_obs (cells)", "n_vars (genes)", "X nnz", "X dtype",
                 "obsm", "layers", "uns keys (sample)"],
        "value": [
            f"{n_obs:,}",
            f"{n_vars:,}",
            f"{meta.get('nnz', 0):,}",
            meta.get("x_dtype", "?"),
            ", ".join(meta["obsm_keys"]) or "—",
            ", ".join(meta["layers"]) or "—",
            ", ".join(meta["uns_keys"][:8]) or "—",
        ],
    })
    sections.append(section("Overview", df_to_html(overview)))

    dtype_rows = []
    for c in sorted(obs_cols):
        col = obs[c]
        n_unique = int(col.nunique(dropna=True))
        sample_vals = col.dropna().astype(str).unique()[:4]
        dtype_rows.append({
            "column": c,
            "dtype": str(col.dtype),
            "n_unique": n_unique,
            "n_missing": int(col.isna().sum()),
            "example_values": ", ".join(sample_vals),
        })
    sections.append(section("obs columns", df_to_html(pd.DataFrame(dtype_rows), max_rows=200)))

    if "skipped_high_cardinality" in meta:
        skipped = pd.DataFrame(meta["skipped_high_cardinality"],
                               columns=["column", "cardinality"])
        sections.append(section("Skipped high-cardinality columns (per-cell)",
                                df_to_html(skipped)))

    candidates_html = "<table><tr><th>role</th><th>guess</th><th>other candidates present</th></tr>"
    for role, picked, pool in [
        ("sample / donor", sample_key, candidate_block(obs_cols, list(dict.fromkeys(SAMPLE_CANDIDATES + DONOR_CANDIDATES)))),
        ("cell_type", celltype_key, candidate_block(obs_cols, CELLTYPE_CANDIDATES)),
        ("batch", batch_key, candidate_block(obs_cols, BATCH_CANDIDATES)),
        ("visit / timepoint", visit_key, candidate_block(obs_cols, VISIT_CANDIDATES)),
        ("sex", sex_key, candidate_block(obs_cols, SEX_CANDIDATES)),
        ("age", age_key, candidate_block(obs_cols, AGE_CANDIDATES)),
        ("condition / outcome", cond_key, candidate_block(obs_cols, CONDITION_CANDIDATES)),
    ]:
        others = [c for c in pool if c != picked]
        candidates_html += (
            f"<tr><td>{escape(role)}</td>"
            f"<td><code>{escape(picked or '—')}</code></td>"
            f"<td class='candidates'>{escape(', '.join(others)) or '—'}</td></tr>"
        )
    candidates_html += "</table>"
    sections.append(section("Suggested keys", candidates_html))

    if donor_key:
        donor_counts = obs[donor_key].value_counts()
        donor_stats = pd.DataFrame({
            "metric": ["n donors", "cells/donor min", "cells/donor median",
                       "cells/donor max", "cells/donor mean"],
            "value": [
                f"{donor_counts.size:,}",
                f"{donor_counts.min():,}",
                f"{int(donor_counts.median()):,}",
                f"{donor_counts.max():,}",
                f"{donor_counts.mean():,.0f}",
            ],
        })
        sections.append(section(f"Donor distribution ({donor_key})",
                                df_to_html(donor_stats) +
                                img_tag(make_donor_cells_hist(donor_counts), "cells per donor")))

    if sample_key and sample_key != donor_key:
        sample_counts = obs[sample_key].value_counts()
        sample_stats = pd.DataFrame({
            "metric": ["n samples", "cells/sample min", "cells/sample median",
                       "cells/sample max"],
            "value": [
                f"{sample_counts.size:,}",
                f"{sample_counts.min():,}",
                f"{int(sample_counts.median()):,}",
                f"{sample_counts.max():,}",
            ],
        })
        if donor_key:
            samples_per_donor = obs.groupby(donor_key, observed=True)[sample_key].nunique()
            sample_stats = pd.concat([sample_stats, pd.DataFrame({
                "metric": ["samples/donor min", "samples/donor median", "samples/donor max"],
                "value": [
                    f"{samples_per_donor.min():,}",
                    f"{int(samples_per_donor.median()):,}",
                    f"{samples_per_donor.max():,}",
                ],
            })], ignore_index=True)
        sections.append(section(f"Sample distribution ({sample_key})",
                                df_to_html(sample_stats) +
                                img_tag(make_donor_cells_hist(sample_counts), "cells per sample")))

    if celltype_key:
        ct_counts = obs[celltype_key].value_counts()
        sections.append(section(f"Cell-type composition ({celltype_key})",
                                df_to_html(ct_counts.head(30).rename("cells").to_frame()) +
                                img_tag(make_celltype_bar(ct_counts), "cell types")))
        if "AIFI_L1" in obs_cols and "AIFI_L2" in obs_cols and "AIFI_L3" in obs_cols:
            granularity = pd.DataFrame({
                "level": ["AIFI_L1", "AIFI_L2", "AIFI_L3"],
                "n_categories": [int(obs[c].nunique()) for c in ["AIFI_L1", "AIFI_L2", "AIFI_L3"]],
            })
            sections.append(section("AIFI annotation granularity",
                                    df_to_html(granularity)))

    if donor_key and visit_key and obs[visit_key].nunique() > 1:
        sections.append(section(f"Donor × {visit_key}",
                                img_tag(make_donor_visit_heatmap(obs, donor_key, visit_key),
                                        "donor visit heatmap")))

    extras_html = ""
    for role, key in [("sex", sex_key), ("age", age_key), ("condition", cond_key), ("visit", visit_key)]:
        if not key:
            continue
        n_unique = obs[key].nunique(dropna=True)
        extras_html += f"<h3>{escape(role)} — <code>{escape(key)}</code></h3>"
        if role == "age" and obs[key].dtype.kind in ("i", "f"):
            extras_html += img_tag(make_age_hist(obs[key]), "age hist")
        elif n_unique > args.max_categorical_values:
            counts = obs[key].value_counts(dropna=False).head(args.max_categorical_values)
            extras_html += df_to_html(counts.rename("cells").to_frame())
            extras_html += f"<p class='candidates'>({n_unique - args.max_categorical_values} more values not shown)</p>"
        else:
            counts = obs[key].value_counts(dropna=False)
            extras_html += df_to_html(counts.rename("cells").to_frame())
            extras_html += img_tag(make_value_count_bar(obs[key], f"{role}: {key}"), f"{role}-bar")
    if extras_html:
        sections.append(section("Other annotations", extras_html))

    batch_cov = [sample_key] if sample_key else []
    if batch_key and batch_key not in batch_cov:
        batch_cov.append(batch_key)
    yaml_hint = f"""<p>Paste into <code>data/{escape(args.cohort)}/dataset_info.yaml</code> after reviewing:</p>
<pre>keys:
  sample: {sample_key or '<unset>'}
  cell_type: {celltype_key or '<unset>'}
  batch: {batch_key or '<unset>'}
  batch_covariates: [{', '.join(batch_cov)}]

metadata:
  samples_metadata_cols: [{', '.join(c for c in [age_key, sex_key, cond_key, visit_key, 'subject.cmv', 'subject.bmi', 'subject.race'] if c and c in obs_cols)}]
  accessible_metadata_columns: [{', '.join(c for c in [age_key, sex_key, 'subject.bmi'] if c and c in obs_cols)}]

evaluate:
  trajectory_variable: {age_key or cond_key or '&lt;unset&gt;'}
  root_sample: &lt;pick the youngest / baseline donor from per-donor stats above&gt;
  inverse_trajectory: false</pre>"""
    sections.append(section("Suggested dataset_info.yaml stanza", yaml_hint))

    html = render_html(args.cohort, sections)
    out_path.write_text(html)
    print(f"wrote {out_path} ({out_path.stat().st_size/1024:.1f} KB)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
