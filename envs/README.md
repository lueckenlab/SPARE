# Auxiliary conda environments

The base/shared environment is defined by `environment.yml` at the repo root
(conda env name `sample_representation_benchmark`). Most components run there.

A few methods have dependency stacks that conflict with the base env and must
live in their own environments. This directory holds their lockfiles.

## `pilot_gm_vae` — the `pilot_gm_vae` represent component

PILOT-GM-VAE (patpy `PILOTGMVAE`, backed by the `pilotgm` fork) **cannot** be
installed into the base env: the fork drags in `pilotpy 2.0.6`, whose hard pins
force `numpy → 1.24`, `pandas → 2.0`, `scipy → 1.11`, `scikit-learn → 1.3`,
`scanpy → 1.9`, breaking squidpy, xarray, lifelines, helical, faiss and `phate`
(which `phemd`/`diffusionemd` need). In a *fresh* env pip instead resolves the
newer `pilotpy 2.0.9`, which is happy with a modern stack (numpy 2.4, pandas 2.3,
scanpy 1.11, torch 2.12).

### Reproduce

```bash
conda create -y -n pilot_gm_vae python=3.11 pip
conda run -n pilot_gm_vae pip install -r envs/pilot_gm_vae-lock.txt
```

`pilot_gm_vae-lock.txt` is a fully-pinned `pip freeze`; every version is fixed so
the resolver cannot fall back to the breaking `pilotpy 2.0.6`. It includes the
exact git commits for `pilotgm` (the fork) and `patpy` (`patpy` is installed
`--no-deps`; its remaining imports are satisfied by the fork's stack).

**Torch / CUDA note:** the fork's default install pulls `torch 2.12.0+cu13`, which
fails on the cluster's CUDA 12.9 GPU driver (`NVIDIA driver too old`). The lock
instead pins `torch 2.6.0+cu124` via `--extra-index-url
https://download.pytorch.org/whl/cu124`. If you reproduce on a node with a newer
driver you can use a later cu12x/cu13 torch build. Run with `--cuda 1` on a GPU
node; `--cuda 0` falls back to CPU (much slower: ~3.5 h for 50 epochs on full
combat even on a V100).

### How it's wired into the pipeline

`src/represent/pilot_gm_vae/` uses the Viash `native` engine, so the Nextflow
process must run with this env active. `diffusionemd` and `phemd` run in the
base env. Splitting per-component envs end-to-end is Phase 0 in `PLAN.md`.
