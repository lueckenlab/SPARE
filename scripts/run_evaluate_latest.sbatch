#!/bin/bash
#SBATCH --job-name=eval_latest
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# Re-aggregate + evaluate using the latest per-method CSVs in
# work_<dataset>/, so we can preview results before sample_representation's
# own aggregate_representations step has fired. Idempotent: re-running
# picks up whatever new CSVs have landed since last invocation.
#
# Usage:
#   sbatch --job-name=eval_<id>_latest --export=DATASET=<id> \
#          scripts/run_evaluate_latest.sbatch
#
# Reads dataset_info.yaml for metadata_path, accessible_metadata_columns,
# cell_type_key, benchmark_schema, root_sample, trajectory_variable,
# inverse_trajectory.

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

if [ -z "${DATASET:-}" ]; then
    echo "DATASET env var required (pass --export=DATASET=combat,...)" >&2
    exit 2
fi

DATASET_INFO="data/${DATASET}/dataset_info.yaml"
if [ ! -f "$DATASET_INFO" ]; then
    echo "Missing $DATASET_INFO" >&2
    exit 2
fi

WORK_DIR="work_${DATASET}"
if [ ! -d "$WORK_DIR" ]; then
    echo "Missing $WORK_DIR — sample_representation hasn't started this dataset yet" >&2
    exit 2
fi

OUTPUT_DIR="output_${DATASET}_evaluate_latest"
STAGE_DIR="${OUTPUT_DIR}/_staged_csvs"
AGG_OUTPUT="${OUTPUT_DIR}/${DATASET}_representations_latest.h5ad"

mkdir -p logs "$OUTPUT_DIR" "$STAGE_DIR"

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

# Stage the newest CSV per basename via a Python helper — symlinks live in
# $STAGE_DIR/ so aggregate_representations sees a clean, dedup'd flat list
# and re-running this sbatch picks up newer mtimes automatically.
python - "$WORK_DIR" "$STAGE_DIR" <<'PY'
import os, sys
from pathlib import Path

work_dir = Path(sys.argv[1])
stage_dir = Path(sys.argv[2])

latest = {}
for csv in work_dir.glob("*/*/*.csv"):
    if not csv.is_file():
        continue
    name = csv.name
    prev = latest.get(name)
    if prev is None or csv.stat().st_mtime > prev.stat().st_mtime:
        latest[name] = csv

# Refresh symlinks
for old in stage_dir.iterdir():
    if old.is_symlink() or old.is_file():
        old.unlink()

for name, src in sorted(latest.items()):
    (stage_dir / name).symlink_to(src.resolve())

print(f"Staged {len(latest)} distinct method CSVs at {stage_dir}")
PY

# Pull metadata / evaluate fields out of dataset_info.yaml. dataset-relative
# paths (files.metadata, evaluate.benchmark_schema) get resolved against the
# data/<id>/ directory.
read -r METADATA_PATH BENCHMARK_SCHEMA CELL_TYPE_KEY ROOT_SAMPLE TRAJECTORY_VARIABLE INVERSE_TRAJECTORY \
    <<<"$(python - "$DATASET_INFO" <<'PY'
import sys, yaml, os
from pathlib import Path

info = yaml.safe_load(open(sys.argv[1]))
yaml_dir = Path(sys.argv[1]).parent

def rel(field):
    val = info.get("files", {}).get(field)
    return str((yaml_dir / val).resolve()) if val else ""

def ev_path(field):
    val = info.get("evaluate", {}).get(field)
    return str((yaml_dir / val).resolve()) if val else ""

metadata_path = rel("metadata")
schema = ev_path("benchmark_schema")
cell_type_key = info.get("keys", {}).get("cell_type") or ""
root_sample = info.get("evaluate", {}).get("root_sample") or ""
traj = info.get("evaluate", {}).get("trajectory_variable") or ""
inverse = "true" if info.get("evaluate", {}).get("inverse_trajectory") else "false"

print(" ".join(f"{x or '__NONE__'}" for x in [metadata_path, schema, cell_type_key, root_sample, traj, inverse]))
PY
)"

if [ "$METADATA_PATH" = "__NONE__" ] || [ ! -f "$METADATA_PATH" ]; then
    echo "metadata_path not resolvable from $DATASET_INFO" >&2
    exit 2
fi
if [ "$BENCHMARK_SCHEMA" = "__NONE__" ] || [ ! -f "$BENCHMARK_SCHEMA" ]; then
    echo "benchmark_schema not resolvable from $DATASET_INFO" >&2
    exit 2
fi

# accessible_metadata_columns gets passed as repeated --accessible_metadata_columns flags
ACCESSIBLE_ARGS=()
while IFS= read -r col; do
    [ -z "$col" ] && continue
    ACCESSIBLE_ARGS+=( --accessible_metadata_columns "$col" )
done < <(python -c "
import sys, yaml
info = yaml.safe_load(open('$DATASET_INFO'))
for c in info.get('metadata', {}).get('accessible_metadata_columns', []) or []:
    print(c)
")

# Build the --input list (multi-arg)
INPUT_ARGS=()
for csv in "$STAGE_DIR"/*.csv; do
    INPUT_ARGS+=( --input "$csv" )
done

echo "=== aggregate_representations ==="
echo "  ${#INPUT_ARGS[@]}/2 method CSVs (--input flags)"
echo "  metadata: $METADATA_PATH"
echo "  cell_type_key: $CELL_TYPE_KEY"
echo "  accessible cols: ${#ACCESSIBLE_ARGS[@]}/2"
echo "  output: $AGG_OUTPUT"

./target/executable/aggregate_representations/aggregate_representations \
    "${INPUT_ARGS[@]}" \
    --output "$AGG_OUTPUT" \
    --metadata_path "$METADATA_PATH" \
    --cell_type_key "$CELL_TYPE_KEY" \
    "${ACCESSIBLE_ARGS[@]}"

echo "=== evaluate ==="
EVAL_ARGS=( \
    --input "$AGG_OUTPUT" \
    --output_dir "$OUTPUT_DIR" \
    --benchmark_schema "$BENCHMARK_SCHEMA" \
    --figure_format png \
)
if [ "$TRAJECTORY_VARIABLE" != "__NONE__" ]; then
    EVAL_ARGS+=( --trajectory_variable "$TRAJECTORY_VARIABLE" )
    EVAL_ARGS+=( --inverse_trajectory "$INVERSE_TRAJECTORY" )
fi
if [ "$ROOT_SAMPLE" != "__NONE__" ]; then
    EVAL_ARGS+=( --root_sample "$ROOT_SAMPLE" )
fi

./target/executable/evaluate/evaluate "${EVAL_ARGS[@]}"

echo "=== done; outputs in $OUTPUT_DIR ==="
ls "$OUTPUT_DIR"
