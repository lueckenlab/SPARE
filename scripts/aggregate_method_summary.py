"""Aggregate benchmark results across datasets.

For each method family and each dataset, pick the *best* representation
(highest `total` score) and plot mean ± std across datasets in the
(relevant, technical) plane — i.e. information retention (X) vs batch
effect removal (Y), matching the existing figure style.

Usage:
    python scripts/aggregate_method_summary.py \
        [--inputs output_*_evaluate_latest/averaged_scores.csv ...] \
        [--metric total] \
        [--output_dir output_method_summary]

Discovers `output_<dataset>_evaluate_latest/averaged_scores.csv` if no
--inputs are given. Outputs:
    output_method_summary/per_method_dataset.csv  -- best rep per
        (method, dataset)
    output_method_summary/aggregated.csv          -- mean / std per method
    output_method_summary/methods_summary.png     -- scatter+errorbar plot
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# Order matters: longer prefixes first so e.g. `grouped_pseudobulk_*` is not
# captured by `pseudobulk_*`, and `pilot_gm_vae_*` is not captured by
# `pilot_*`.
PREFIX_RULES = [
    ("grouped_pseudobulk", "CT pseudobulk"),
    ("cell_group_composition", "Cell type composition"),
    ("pilot_gm_vae", "PILOT-GM-VAE"),
    ("sampleclr", "SampleCLR"),
    ("scpoli_sample_embeddings", "scPoli sample embeddings"),
    ("scpoli", "scPoli"),
    ("gloscope", "GloScope"),
    ("mofa", "MOFA"),
    ("pilot", "PILOT"),
    ("pseudobulk", "Pseudobulk"),
    ("phemd", "PhEMD"),
    ("diffusionemd", "DiffusionEMD"),
    ("mrvi", "MrVI"),
    ("random_vector", "Random vector"),
    ("ehrapy", "Ehrapy"),
]

# Handcrafted palette. The first nine entries (GloScope through Random
# vector) match the colours from the figure in the writeup; the remainder
# extend that palette to the methods added since.
METHOD_COLORS = {
    "GloScope": "#D26995",                 # magenta-pink
    "MOFA": "#E5C239",                     # gold/yellow
    "PILOT": "#E48F2A",                    # orange
    "Cell type composition": "#3FA585",    # teal/green
    "Pseudobulk": "#83B0CE",               # light blue
    "CT pseudobulk": "#3F87B5",            # blue
    "scPoli sample embeddings": "#D4642D", # red-orange
    "Ehrapy": "#000000",                   # black
    "Random vector": "#808080",            # gray
    # New methods (extension of the writeup palette, distinct hues)
    "SampleCLR": "#7B3FB5",                # purple
    "MrVI": "#8C564B",                     # brown
    "PILOT-GM-VAE": "#7F4310",             # dark amber (PILOT family, darker)
    "scPoli": "#C13F25",                   # deeper red (scPoli family, darker)
    "PhEMD": "#9C8FA0",                    # muted plum
    "DiffusionEMD": "#5FB57A",             # leaf green
    "Other": "#444444",                    # dark gray catch-all
}


def classify_representation(name: str) -> str:
    for prefix, family in PREFIX_RULES:
        if name == prefix or name.startswith(prefix + "_") or name.startswith(prefix):
            if name == prefix or name.startswith(prefix + "_"):
                return family
    return "Other"


def best_per_method(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Pick the highest-`metric` representation within each method family."""
    df = df.copy()
    df["method"] = df.index.map(classify_representation)
    df["representation"] = df.index
    df = df.sort_values(metric, ascending=False)
    return df.drop_duplicates("method", keep="first").set_index("method")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs",
        nargs="*",
        help="averaged_scores.csv files. Defaults to "
        "output_*_evaluate_latest/averaged_scores.csv.",
    )
    parser.add_argument(
        "--metric",
        default="total",
        help="Score used to pick the best representation within each "
        "method family per dataset (default: total).",
    )
    parser.add_argument(
        "--output_dir",
        default="output_method_summary",
        help="Directory for the aggregated CSV + figure.",
    )
    parser.add_argument(
        "--figure_format",
        default="png",
        choices=("png", "pdf", "svg"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inputs = args.inputs or sorted(
        Path(".").glob("output_*_evaluate_latest/averaged_scores.csv")
    )
    if not inputs:
        print("No averaged_scores.csv files found", file=sys.stderr)
        sys.exit(2)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    per_method_dataset = []
    for path in inputs:
        path = Path(path)
        match = re.match(
            r"output_(?P<id>.+?)_evaluate(?:_latest)?", path.parent.name
        )
        dataset = match.group("id") if match else path.parent.name
        df = pd.read_csv(path, index_col=0)
        if args.metric not in df.columns:
            print(
                f"  skip {dataset}: missing column {args.metric!r}",
                file=sys.stderr,
            )
            continue
        best = best_per_method(df, args.metric)
        best["dataset"] = dataset
        per_method_dataset.append(best.reset_index())
        print(
            f"  {dataset}: {len(best)} methods present "
            f"(top representation: {best.iloc[0]['representation']})"
        )

    if not per_method_dataset:
        print("No usable inputs", file=sys.stderr)
        sys.exit(2)

    combined = pd.concat(per_method_dataset, ignore_index=True)

    metrics = [
        c for c in ("relevant", "technical", "trajectory", "contextual", "total")
        if c in combined.columns
    ]

    # Gain-over-pseudobulk per dataset, per metric.
    # gain = (score - pb) / (1 - pb), with pb = max Pseudobulk score on that
    # metric in that dataset. 0 means parity with pseudobulk, 1 means perfect,
    # negative means worse. NaN if no Pseudobulk row landed for the dataset
    # or pb == 1 (degenerate).
    pb_rows = combined[combined["method"] == "Pseudobulk"]
    if pb_rows.empty:
        print(
            "WARNING: no Pseudobulk rows found — gain columns will be NaN",
            file=sys.stderr,
        )
    for metric in metrics:
        gain_col = f"{metric}_gain"
        combined[gain_col] = float("nan")
        for dataset, group in combined.groupby("dataset"):
            pb_score = pb_rows[pb_rows["dataset"] == dataset][metric].max()
            if pd.isna(pb_score) or pb_score >= 1.0:
                continue
            combined.loc[group.index, gain_col] = (
                group[metric] - pb_score
            ) / (1.0 - pb_score)

    combined_path = out_dir / "per_method_dataset.csv"
    combined.to_csv(combined_path, index=False)
    print(f"Wrote {combined_path}")

    agg_cols = metrics + [f"{m}_gain" for m in metrics]
    aggregated = (
        combined.groupby("method")[agg_cols]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    aggregated.columns = [
        "_".join(filter(None, col)) if isinstance(col, tuple) else col
        for col in aggregated.columns
    ]
    aggregated = aggregated.sort_values("total_mean", ascending=False)
    aggregated_path = out_dir / "aggregated.csv"
    aggregated.to_csv(aggregated_path, index=False)
    print(f"Wrote {aggregated_path}")

    method_order = (
        aggregated.sort_values("total_mean", ascending=False)["method"].tolist()
    )
    color_for = lambda m: METHOD_COLORS.get(m, METHOD_COLORS["Other"])

    def _scatter_plot(x_col: str, y_col: str, x_label: str, y_label: str,
                       out_name: str, *, gain_axes: bool = False) -> None:
        if not {x_col, y_col}.issubset(combined.columns):
            print(
                f"averaged_scores lacks {x_col!r} or {y_col!r}; skipping "
                f"{out_name}.",
                file=sys.stderr,
            )
            return
        # Pseudobulk is the reference line in gain plots; don't draw it (it
        # sits at (0,0) by construction). On absolute-score plots we still
        # show it like any other method.
        plot_methods = [
            m for m in method_order
            if not (gain_axes and m == "Pseudobulk")
        ]
        fig, ax = plt.subplots(figsize=(7, 7))
        for method in plot_methods:
            group = combined[combined["method"] == method].dropna(
                subset=[x_col, y_col]
            )
            if len(group) == 0:
                continue
            x_mean = group[x_col].mean()
            y_mean = group[y_col].mean()
            x_std = group[x_col].std(ddof=1) if len(group) > 1 else 0.0
            y_std = group[y_col].std(ddof=1) if len(group) > 1 else 0.0
            ax.errorbar(
                x_mean, y_mean, xerr=x_std, yerr=y_std,
                fmt="o",
                color=color_for(method),
                ecolor=color_for(method),
                elinewidth=1.2,
                capsize=3,
                markersize=8,
                label=f"{method} (n={len(group)})",
            )
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        if gain_axes:
            # Pseudobulk reference at (0, 0); show a clear cross there.
            ax.axhline(0, color="lightgray", lw=1, zorder=0)
            ax.axvline(0, color="lightgray", lw=1, zorder=0)
            ax.set_xlim(-1, 1)
            ax.set_ylim(-1, 1)
        else:
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
        ax.grid(False)
        sns.despine(ax=ax)
        ax.legend(
            title="Method", loc="center left", bbox_to_anchor=(1.02, 0.5),
            frameon=False,
        )
        fig.tight_layout()
        path = out_dir / f"{out_name}.{args.figure_format}"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {path}")

    _scatter_plot(
        "relevant", "technical",
        "Information retention (relevant)",
        "Batch effect removal (technical)",
        "methods_summary",
    )
    _scatter_plot(
        "relevant", "trajectory",
        "Information retention (relevant)",
        "Trajectory preservation",
        "methods_relevant_vs_trajectory",
    )
    _scatter_plot(
        "relevant_gain", "technical_gain",
        "Relevant gain vs Pseudobulk",
        "Technical gain vs Pseudobulk",
        "methods_gain_summary",
        gain_axes=True,
    )
    _scatter_plot(
        "relevant_gain", "trajectory_gain",
        "Relevant gain vs Pseudobulk",
        "Trajectory gain vs Pseudobulk",
        "methods_gain_relevant_vs_trajectory",
        gain_axes=True,
    )


if __name__ == "__main__":
    main()
