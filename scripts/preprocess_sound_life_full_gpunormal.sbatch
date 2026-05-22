#!/bin/bash
#SBATCH --job-name=pp_sl_norm
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --mem=1300G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/pp_sl_norm_%j.out
#SBATCH --error=logs/pp_sl_norm_%j.err
#
# Sole production preprocess of the sound_life cohort (13.79M cells) on
# gpu_normal (24h cap, NON-preemptible). Moved here from gpu_preemptible
# after a preemptible run was bumped off its node mid-scANVI: with scANVI
# capped at 3 epochs (--n_scanvi_epochs 3) the full run is only ~5-6h,
# well under the 24h cap, so the preemptible queue's restart-from-scratch
# risk is no longer worth it. Carries the NaN-stability fix
# (--gradient_clip_val 1.0; see preprocess_sound_life_full.sbatch header
# for the full failure history). Writes the canonical output paths.

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
