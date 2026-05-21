#!/bin/bash
#SBATCH --job-name=pp_mi
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --nice=10000
#SBATCH --output=logs/pp_mi_%j.out
#SBATCH --error=logs/pp_mi_%j.err
#
# Preprocess the Kuppe 2022 myocardial infarction dataset (191,795 cells,
# 29 samples, 20 donors). Small relative to the AIFI cohorts, so 128G +
# one GPU on gpu_normal (24h cap) is ample; scVI auto-scales to ~40
# epochs at this cell count and scANVI is fixed at 20.
#
# Covariate names follow the recent-dataset convention (real obs column
# names): scVI/scANVI run once per batch covariate, producing
# X_scVI_sample / X_scVI_donor_id and the matching X_scANVI_* embeddings,
# plus X_pca and X_scpoli.
#
# Runs: HVG -> full-library normalize -> PCA -> scVI/scANVI per batch
# covariate (sample, donor_id) -> scPoli -> QC + metadata.

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs

SECONDS=0
./target/executable/preprocess/preprocess \
    --input  data/myocardial_infarction/myocardial_infarction_cleaned.h5ad \
    --output data/myocardial_infarction/myocardial_infarction_processed.h5ad \
    --output_metadata data/myocardial_infarction/myocardial_infarction_metadata.csv \
    --sample_key "sample" \
    --cell_type_key "cell_type_original" \
    --batch_key "donor_id" \
    --batch_covariates "sample,donor_id" \
    --samples_metadata_cols "donor_id,patient_region_id,patient_group,major_labl,disease,sex,development_stage,age,disease_progression" \
    --sample_size_threshold 100 \
    --output_compression gzip
echo "Wallclock: ${SECONDS}s"
