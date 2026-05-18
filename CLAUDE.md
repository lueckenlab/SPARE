# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A benchmark of single-cell sample/patient representation methods, built as a **Viash + Nextflow** pipeline. Components (`src/<stage>/<name>/`) compile to executables and Nextflow modules under `target/`. Workflows (`src/workflows/<name>/`) wire components into a DAG. The flow is: `download_data → clean_data → preprocess → represent (×N methods) → aggregate_representations → evaluate`.

Methods compared: `pseudobulk`, `grouped_pseudobulk`, `cell_group_composition`, `random` (baseline), `scpoli`, `gloscope`, `mofa`, `mrvi`, `pilot`. Each writes a sample×sample distance matrix or sample embedding as CSV; `aggregate_representations` collects them into a single AnnData where each representation is `.obsm[<method>_<experiment>]` and a metadata-only ehrapy representation is built from `--accessible_metadata_columns`. `evaluate` reads that AnnData and computes signal retention, batch removal, trajectory preservation, replicate robustness, scalability.

Datasets ingested via `src/clean_data/<dataset>/`: `combat`, `stephenson`, `hlca`, `onek1k`, `copd`, `ticatlas`. Each has dataset-specific cleaning (annotation fixes, metadata reformatting) before the shared `preprocess` (normalize → log → HVG → PCA → Harmony → scVI → scANVI → QC → metadata extraction).

Note: `src/preprocess/` is still driven by its own `params.yaml` (legacy interface); it has not been migrated to `dataset_info.yaml` yet. Other workflows (download_datasets, clean_data, sample_representation) use the new pattern.

## Build & run

The Conda env `spare` (defined by `environment.yml` at the repo root) is required. Activate it before any viash/nextflow command. patpy is pinned to commit `44903dc` — the latest tagged release (0.15.2) is behind the working version (0.16.5), and we depend on `GloScope_py` + `patpy.tl.supervised` which only exist past the last tag.

Note: `nextflow.config` currently sets `LD_LIBRARY_PATH` to a hard-coded user path. This will be removed; until then, anyone reproducing the pipeline needs to edit that line or override the env var.

```bash
# Build everything in src/ → target/
viash ns build --setup cb

# Build a single component
viash build src/represent/scpoli/config.vsh.yaml -p nextflow -o target/nextflow/scpoli

# Run one component standalone (uses the executable runner)
./target/executable/preprocess/preprocess --input data/combat/combat_cleaned.h5ad --output out.h5ad ...

# Run a sub-workflow via Nextflow (from repo root)
nextflow run target/nextflow/sample_representation/main.nf \
  -profile standard \
  -params-file src/workflows/sample_representation/datasets.yaml \
  --method_params src/workflows/sample_representation/experiments.yaml \
  --publish_dir output/

# Run tests for one component (writes outputs under data/)
viash test src/preprocess/config.vsh.yaml
```

The cluster runs Slurm. `nextflow.config` defines profiles that submit `gpu`, `highcpu`, `singlecpu` labelled processes to the appropriate queues (`gpu_p`, `cpu_p`). `process.errorStrategy` retries on OOM-related exit codes (137–140) and ignores other failures, so a single bad method doesn't abort the matrix. `src/workflows/utils/labels.config` defines memory/CPU labels (`lowmem`/`midmem`/`highmem`/`veryhighmem` and `singlecpu`/`lowcpu`/`midcpu`/`highcpu`) consumed by per-component `directives.label` in `config.vsh.yaml`.

## Anatomy of a component

Every component is a directory with:
- `config.vsh.yaml` — Viash schema: arguments, resources, engines (`native`), runners (`executable` + `nextflow`), Nextflow `directives.label` for cluster resources.
- `script.py` — the actual logic. Top of file has a `## VIASH START` / `## VIASH END` block defining `par = {...}` so the script is runnable standalone in a notebook; Viash replaces this block with real CLI args at build time.
- Optional `test.py` referenced as `test_resources`.

When adding a new representation method:
1. Create `src/represent/<name>/{config.vsh.yaml, script.py}` — write CSV of (sample × sample distances) or (sample × dims) indexed by sample ID.
2. Add `<name>` to the `methods` list in `src/workflows/sample_representation/main.nf`.
3. Add `<name>` as a dependency in `src/workflows/sample_representation/config.vsh.yaml`.
4. Add an entry under `<name>:` in `experiments.yaml` mapping experiment-name → params dict (these become CLI args for the script).
5. `viash ns build --setup cb` to rebuild `target/`.

## Per-dataset config

All dataset-specific parameters live in **`data/<id>/dataset_info.yaml`** (canonical schema: `src/workflows/utils/dataset_info_example.yaml`). The file defines `files.*` paths, `keys.*` (sample/cell_type/batch column names), `metadata.{samples_metadata_cols, accessible_metadata_columns}`, plus `preprocess.*` and `evaluate.*` settings. The dataset YAMLs are git-tracked via a `.gitignore` exception.

Each workflow takes `--dataset_info path/to/dataset_info.yaml` as its primary input. The Groovy helper `src/workflows/utils/dataset_info.nf::expandDatasetInfo` projects the YAML into a flat state map keyed to Viash component arg names. Individual CLI args (`--input`, `--sample_key`, etc.) are still accepted and **win when explicitly set** — useful for ad-hoc overrides.

Per-workflow `datasets.yaml` registry files list which datasets to iterate over:
```yaml
- id: combat
  dataset_info: ../../../data/combat/dataset_info.yaml
```

## How `sample_representation` expands the experiment matrix

`src/workflows/sample_representation/main.nf` reads `--method_params` (an `experiments.yaml`), then for each dataset in `datasets.yaml` it flatMaps over `methods × experiments[method]`, creating a run id `${dataset}.${method}.${experiment}`. `runEach` filters so each run executes only the matching component. After all methods finish, results are grouped back by original dataset id and fed to `aggregate_representations`. So `experiments.yaml` controls *which method+layer combinations are evaluated*; `datasets.yaml` controls *which datasets they run on*.

Test runs: use `experiments_test.yaml` (tiny epochs, single layers) with the same `datasets.yaml`. Subsampled variants are selected by pointing `--dataset_info` at a YAML whose `files.processed` references a subsampled h5ad.

## Notable conventions in scripts

- All scripts use `random.seed(42)` and `np.random.seed(42)` at the top — preserve this.
- `preprocess` always stores raw counts in `adata.layers["X_raw_counts"]` and log-normalized in `adata.layers["X_log_norm_counts"]`. Downstream methods reference layers by these names (see `experiments.yaml`'s `layer:` keys: `X_raw_counts`, `X_log_norm_counts`, `X_pca`, `X_scpoli`, `X_scVI_sample`, `X_scVI_batch`, `X_scANVI_sample`, `X_scANVI_batch`).
- `aggregate_representations/script.py` clips negative values to 0 (gloscope produces them) and skips any representation with >10% NaNs.
- `accessible_metadata_columns` (formerly `obs_columns` — renamed in commit `360dac3`) defines the *measurable clinical features* used to build the ehrapy baseline representation. Everything not in that list goes to `.obs` for evaluation.
- Sample identifiers must be unique (commit `bd50c4d` enforces this via `obs_names_make_unique()`).

## Plotting

Follow the global rules in `~/.claude/CLAUDE.md`: no grid lines, `sns.despine()`, no frame/ticks on UMAPs (`frameon=False`), normalized metrics get `ylim=(0, 1)`, signed metrics get symmetric range with `RdBu_r`.
