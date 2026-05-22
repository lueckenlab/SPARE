#!/bin/bash
#SBATCH --job-name=pp_sl_full
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_preemptible
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --mem=1300G
#SBATCH --time=3-00:00:00
#SBATCH --nice=10000
#SBATCH --output=logs/pp_sl_full_%j.out
#SBATCH --error=logs/pp_sl_full_%j.err
#
# Full production preprocess of the sound_life cohort (13.79M cells,
# no subsampling). The cleaned file is 259GB. Observed peak RSS across
# two earlier runs was 774GB (stable through scANVI), so 1.3TB gives
# ~1.7x headroom and fits on far more partially-used 2TB nodes than the
# original 1.8TB request, which sat at the edge and starved scheduling.
#
# STABILITY: scVI/scANVI diverge to NaN (encoder emits nan loc in
# Normal()) at random epochs on this 13.8M-cell cohort -- a gradient
# explosion. Lowering epochs did NOT fix it (20-epoch runs died @8/@17,
# 5-epoch @4, 3-epoch @3 on the 2nd covariate). The real fix is
# --gradient_clip_val 1.0, which clips the grad norm on every scVI+scANVI
# train call (off by default, so other datasets are unaffected).
# --n_scanvi_epochs kept at 3 here to shorten the run. NOTE:
# gpu_preemptible jobs can be killed if a priority job needs the node,
# and the preprocess script has no checkpoint/resume -- if preempted,
# the run must restart from scratch.
#
# Runs: HVG -> full-library normalize -> PCA -> scVI/scANVI per batch
# covariate (sample.sampleKitGuid, batch_id) -> scPoli -> QC + metadata.

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
    --output_compression gzip
echo "Wallclock: ${SECONDS}s"
