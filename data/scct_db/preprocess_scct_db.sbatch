#!/bin/bash
#SBATCH --job-name=pp_scct
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/pp_%x_%j.out
#SBATCH --error=logs/pp_%x_%j.err
#
# Preprocess one CD45 compartment of the scCT-DB GSE246613 cohort.
# Usage: sbatch -J pp_scct_db_cd45pos data/scct_db/preprocess_scct_db.sbatch scct_db_cd45pos
#
# No batch covariate: batch_covariates=[sample_id] -> only X_*_sample
# embeddings (+ X_pca, X_scpoli). batch_key=source is a constant column so
# HVG is global (no real batch). `dataset` is never a covariate (it encodes
# the response label).

set -euo pipefail
DS="${1:?usage: preprocess_scct_db.sbatch <scct_db_cd45pos|scct_db_cd45neg>}"
cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark
source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark
mkdir -p logs

SECONDS=0
./target/executable/preprocess/preprocess \
    --input  "data/${DS}/${DS}_cleaned.h5ad" \
    --output "data/${DS}/${DS}_processed.h5ad" \
    --output_metadata "data/${DS}/${DS}_metadata.csv" \
    --sample_key "sample_id" \
    --cell_type_key "cell_type" \
    --batch_key "source" \
    --batch_covariates "sample_id" \
    --samples_metadata_cols "response,timepoint,sort,dataset,patient,source" \
    --sample_size_threshold 100 \
    --n_scanvi_epochs 20 \
    --gradient_clip_val 1.0 \
    --output_compression gzip
echo "Wallclock: ${SECONDS}s  ($DS)"
