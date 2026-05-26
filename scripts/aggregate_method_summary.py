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
    combined_path = out_dir / "per_method_dataset.csv"
    combined.to_csv(combined_path, index=False)
    print(f"Wrote {combined_path}")

    metrics = [
        c for c in ("relevant", "technical", "trajectory", "contextual", "total")
        if c in combined.columns
    ]
    aggregated = (
        combined.groupby("method")[metrics]
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

    if not {"relevant", "technical"}.issubset(combined.columns):
        print(
            "averaged_scores files lack 'relevant' or 'technical' columns; "
            "skipping plot.",
            file=sys.stderr,
        )
        return

    palette = sns.color_palette("husl", n_colors=combined["method"].nunique())
    method_order = (
        aggregated.sort_values("total_mean", ascending=False)["method"].tolist()
    )
    color_map = dict(zip(method_order, palette))

    fig, ax = plt.subplots(figsize=(7, 7))
    for method in method_order:
        group = combined[combined["method"] == method]
        if len(group) == 0:
            continue
        rel_mean = group["relevant"].mean()
        rel_std = group["relevant"].std(ddof=1) if len(group) > 1 else 0.0
        tech_mean = group["technical"].mean()
        tech_std = group["technical"].std(ddof=1) if len(group) > 1 else 0.0
        ax.errorbar(
            rel_mean,
            tech_mean,
            xerr=rel_std,
            yerr=tech_std,
            fmt="o",
            color=color_map[method],
            ecolor=color_map[method],
            elinewidth=1.2,
            capsize=3,
            markersize=8,
            label=f"{method} (n={len(group)})",
        )

    ax.set_xlabel("Information retention (relevant)")
    ax.set_ylabel("Batch effect removal (technical)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(False)
    sns.despine(ax=ax)
    ax.legend(
        title="Method", loc="center left", bbox_to_anchor=(1.02, 0.5),
        frameon=False,
    )

    fig.tight_layout()
    figure_path = out_dir / f"methods_summary.{args.figure_format}"
    fig.savefig(figure_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {figure_path}")


if __name__ == "__main__":
    main()
