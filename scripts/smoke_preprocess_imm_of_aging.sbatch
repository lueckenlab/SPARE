#!/bin/bash
#SBATCH --job-name=pp_ioa_smoke
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=8:00:00
#SBATCH --nice=10000
#SBATCH --output=logs/pp_ioa_smoke_%j.out
#SBATCH --error=logs/pp_ioa_smoke_%j.err
#
# Smoke test of src/preprocess on a stratified subsample of the
# imm_of_aging cleaned data. Builds the subsample first (backed-mode,
# cheap) if missing, then runs the full preprocess (HVG -> normalize ->
# PCA -> scVI/scANVI per batch covariate -> scPoli -> metadata).
#
# Purpose: confirm the AIFI keys flow through and the HVG/normalize
# refactor works end-to-end on real data. NOT a production run.

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs

SMOKE=data/imm_of_aging/imm_of_aging_smoke.h5ad
if [[ ! -f "${SMOKE}" ]]; then
    echo "=== Building subsample ${SMOKE} ==="
    python scripts/subsample_h5ad.py \
        --input  data/imm_of_aging/imm_of_aging_cleaned.h5ad \
        --output "${SMOKE}" \
        --sample-key sample.sampleKitGuid \
        --per-sample 1300
fi

echo "=== Running preprocess ==="
SECONDS=0
./target/executable/preprocess/preprocess \
    --input  "${SMOKE}" \
    --output data/imm_of_aging/imm_of_aging_smoke_processed.h5ad \
    --output_metadata data/imm_of_aging/imm_of_aging_smoke_metadata.csv \
    --sample_key "sample.sampleKitGuid" \
    --cell_type_key "AIFI_L2" \
    --batch_key "batch_id" \
    --batch_covariates "sample.sampleKitGuid,batch_id" \
    --samples_metadata_cols "subject.ageAtFirstDraw,subject.biologicalSex,subject.cmv,subject.race,subject.ethnicity,sample.visitName,cohort.cohortGuid,batch_id,pool_id,chip_id" \
    --sample_size_threshold 100 \
    --output_compression gzip
echo "Wallclock: ${SECONDS}s"
