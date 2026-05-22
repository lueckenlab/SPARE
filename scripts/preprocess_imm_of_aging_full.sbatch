#!/bin/bash
#SBATCH --job-name=pp_ioa_full
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=400G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/pp_ioa_full_%j.out
#SBATCH --error=logs/pp_ioa_full_%j.err
#
# Full production preprocess of the imm_of_aging cohort (3.76M cells,
# no subsampling). 400G + GPU lands on a supergpu (2TB) node so the
# 111GB cleaned file loads and scVI/scANVI/scPoli run end-to-end.
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
    --input  data/imm_of_aging/imm_of_aging_cleaned.h5ad \
    --output data/imm_of_aging/imm_of_aging_processed.h5ad \
    --output_metadata data/imm_of_aging/imm_of_aging_metadata.csv \
    --sample_key "sample.sampleKitGuid" \
    --cell_type_key "AIFI_L2" \
    --batch_key "batch_id" \
    --batch_covariates "sample.sampleKitGuid,batch_id" \
    --samples_metadata_cols "subject.ageAtFirstDraw,subject.biologicalSex,subject.cmv,subject.race,subject.ethnicity,sample.visitName,cohort.cohortGuid,batch_id,pool_id,chip_id" \
    --sample_size_threshold 100 \
    --output_compression gzip
echo "Wallclock: ${SECONDS}s"
