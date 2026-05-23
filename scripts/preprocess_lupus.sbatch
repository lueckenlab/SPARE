#!/bin/bash
#SBATCH --job-name=pp_lupus
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=07:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/pp_lupus_%j.out
#SBATCH --error=logs/pp_lupus_%j.err
#
# Step 2/2 of the lupus pipeline: preprocess (scVI/scANVI/scPoli + QC).
# lupus is 1.2M cells (~combat-scale), so we use preprocess defaults
# (n_scanvi_epochs=20, n_scpoli_epochs=50, no gradient clipping) — none
# of the sound_life-scale knobs are needed at this size.

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs

SECONDS=0
./target/executable/preprocess/preprocess \
    --input  data/lupus/lupus_cleaned.h5ad \
    --output data/lupus/lupus_processed.h5ad \
    --output_metadata data/lupus/lupus_metadata.csv \
    --sample_key "ind_cov_batch_cov" \
    --cell_type_key "cg_cov" \
    --batch_key "Processing_Cohort" \
    --batch_covariates "ind_cov_batch_cov,Processing_Cohort" \
    --samples_metadata_cols "Processing_Cohort,Age,Sex,pop_cov" \
    --sample_size_threshold 100 \
    --output_compression gzip
echo "Wallclock: ${SECONDS}s"
