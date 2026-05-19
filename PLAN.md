# SPARE development plan

Tracks the multi-phase refactor for the patient/sample representation benchmark.
Update as items land; commit when ticked.

Open GitHub issues:
- [#1](https://github.com/lueckenlab/pat_rep_benchmark/issues/1) improve HLCA cleaning
- [#2](https://github.com/lueckenlab/pat_rep_benchmark/issues/2) improve metadata processing
- [#7](https://github.com/lueckenlab/pat_rep_benchmark/issues/7) rename replicate-robustness metric → FOSCTTM
- [#9](https://github.com/lueckenlab/pat_rep_benchmark/issues/9) supervised baseline
- [#10](https://github.com/lueckenlab/pat_rep_benchmark/issues/10) add Human Immune Cell Atlas

---

## Phase 0 — Prerequisites

- [ ] Review and merge open PRs (dataset PR + foundation-model PR).
- [x] Pin patpy version. Pinned to commit `44903dc` in `environment.yml`.
- [ ] Capture conda envs. Currently a single `spare` env. Split into:
  - `spare-cpu` — scanpy/scvi-tools/patpy core (preprocess, classical methods, aggregate, evaluate).
  - `spare-gpu` — torch + scGPT deps (mrvi, scpoli, preprocess_scgpt).
  - `spare-uce` — UCE has conflicting deps, isolate.
  - `spare-supervised` — `patpy.tl.supervised` deps (mixmil, pulsar). Reuse `spare-gpu` if compatible.
- [ ] Drop hard-coded `LD_LIBRARY_PATH` from `nextflow.config` (line 17 points at one user's home).

## Phase 1 — Per-dataset YAML consolidation (issue #2)

- [x] `data/<id>/dataset_info.yaml` schema defined and adopted for `download_datasets`, `clean_data`, `sample_representation`. Canonical example: `src/workflows/utils/dataset_info_example.yaml`. Groovy helper at `src/workflows/utils/dataset_info.nf::expandDatasetInfo`.
- [ ] Migrate `src/preprocess/` off its legacy `params.yaml` onto `dataset_info.yaml`.

## Phase 2 — Patpy alignment refactor

- [ ] Replace inline scanpy in `src/preprocess/script.py` with `patpy.pp.basic` primitives where they exist.
- [ ] Audit each `src/represent/*/script.py` against the pinned patpy API; upstream fixes go to patpy, not local wrappers.
- [ ] Issue #7 — rename evaluate metric → FOSCTTM.
- [ ] Issue #1 — incorporate HLCA cleaning improvements.

## Phase 3 — Foundation-model preprocessing

FM components produce a single `adata.obsm["X_<fm>"]` layer; only pseudobulk consumes it initially.

- [ ] `src/preprocess_scgpt/` — `spare-gpu`, GPU label.
- [ ] `src/preprocess_uce/` — `spare-uce`, GPU label.
- [ ] `experiments.yaml` additions:
  ```yaml
  pseudobulk:
    scgpt_euclidean: { layer: X_scgpt, distance: euclidean }
    scgpt_cosine:    { layer: X_scgpt, distance: cosine }
    uce_euclidean:   { layer: X_uce,   distance: euclidean }
    uce_cosine:      { layer: X_uce,   distance: cosine }
  ```
- [ ] Decide: write extra `obsm` layers into the existing `processed.h5ad`, or sidecar files.

## Phase 4 — Supervised baseline (issue #9)

- [ ] `src/represent/supervised/` — one component per method or one parametrized using patpy's `MixMIL`, `PULSAR`, `PaSCient` via `SupervisedSampleMethod`.
- [ ] Contract: writes sample×dims embedding CSV like other methods, plus `--label_keys`, `--tasks` (classification/regression/ranking).
- [ ] Decide train/test split location (in component vs. evaluate).
- [ ] Add `supervised:` block in `experiments.yaml`.
- [ ] Evaluate metrics for supervised differ from unsupervised — flag for evaluate refactor.

## Phase 5 — New datasets

- [x] TICAtlas — added (`src/clean_data/ticatlas/`, `data/ticatlas/dataset_info.yaml`, full launcher in `scripts/run_ticatlas_full.sbatch`).
- [ ] SEA-AD (Seattle Alzheimer's, brain) — public on CELLxGENE.
- [ ] HICA (issue #10) — Human Immune Cell Atlas.
- [ ] Immunobiology of aging (`aifi_data/imm_of_aging`) — counts in `.raw.X`, age `89+→89`. Reuse `patpy/scripts/preprocess_aifi.py`.
- [ ] Vaccine prediction (`aifi_data/sound_life`?) — confirm cohort.
- [ ] Open-PR dataset — depends on PR review.

## Phase 6 — End-to-end on all datasets

- [x] `process.errorStrategy = { task.exitStatus in 137..140 ? 'retry' : 'ignore' }` — retries on OOM, ignores other failures so the matrix keeps going.
- [ ] Status-collector component at end of `sample_representation` — writes `summary.csv` of (dataset, method, experiment, status, runtime, peak memory).
- [ ] Per-dataset memory tuning: HLCA, SEA-AD likely need bumped `veryhighmem`. Brain atlases may need a subsample profile.
- [ ] `runIf: { state.cleaned_file && state.cleaned_file.exists() }` on download components so the pipeline auto-skips when cleaned data is present.

## Phase 7 — Docs & repro

- [ ] Update `README.md` with "How to run" beyond the mermaid.
- [ ] `data/README.md` summarizing each dataset.
- [ ] Update `CLAUDE.md` for new YAML schema and FM/supervised components.
- [ ] Wire `evaluate` into a workflow (currently a standalone component; `src/workflows/run_pipeline/main.nf` has TODOs).

## Cross-cutting: patpy-direction policy

While refactoring, if a SPARE script has a helper that's genuinely general (cleaning recipes, QC heuristics, metric implementations), move it to patpy and open a PR there. Conversely, if patpy has an API that's awkward for SPARE, fix it upstream rather than wrapping it locally.

## Loose ends (small, not phased)

- [ ] `mofa.log_counts` fails on TICAtlas — references `X_log_norm_counts` but the ticatlas processed h5ad only has `X_raw_counts` and `shifted_log_counts`. Currently swallowed by `errorStrategy=ignore`. Fix in `experiments.yaml` or re-run preprocess.
