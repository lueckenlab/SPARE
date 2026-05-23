#!/bin/bash
#SBATCH --job-name=pp_sl_sp5
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --mem=1300G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/pp_sl_sp5_%j.out
#SBATCH --error=logs/pp_sl_sp5_%j.err
#
# Restart preprocess of sound_life with scPoli capped at 5 epochs.
# The first run (job 36661881) hit a wall: at the observed ~15 min/percent
# scPoli rate, projected completion overran the 24h gpu_normal cap by ~6h,
# meaning the in-memory state (no checkpoints in preprocess) would be lost.
# This run uses the newly-added --n_scpoli_epochs knob (default 50) to cap
# scPoli at 5 — should bring scPoli wall time from ~12h to ~1h, leaving
# plenty of headroom even if scVI/scANVI take their usual ~3h together.
#
# We deliberately do NOT cancel the old run: if it somehow finishes before
# this one, we cancel this one and keep the 50-epoch result. Otherwise it
# dies at its 24h cap and this run wins.

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs

SECONDS=0
./target/executable/preprocess/preprocess \
    --input  data/sound_life/sound_life_cleaned.h5ad \
    --output data/sound_life/sound_life_processed.h5ad \
    --output_metadata data/sound_life/sound_life_metadata.csv \
    --sample_key "sample.sampleKitGuid" \
    --cell_type_key "AIFI_L2" \
    --batch_key "batch_id" \
    --batch_covariates "sample.sampleKitGuid,batch_id" \
    --samples_metadata_cols "subject.ageAtFirstDraw,subject.ageGroup,subject.biologicalSex,subject.bmi,subject.cmv,subject.race,subject.ethnicity,sample.visitName,sample.drawYear,cohort.cohortGuid,batch_id,pool_id,chip_id" \
    --sample_size_threshold 100 \
    --n_scanvi_epochs 3 \
    --gradient_clip_val 1.0 \
    --n_scpoli_epochs 5 \
    --output_compression gzip
echo "Wallclock: ${SECONDS}s"
