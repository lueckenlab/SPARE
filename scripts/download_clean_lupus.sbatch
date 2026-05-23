#!/bin/bash
#SBATCH --job-name=dl_clean_lupus
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --nice=10000
#SBATCH --output=logs/dl_clean_lupus_%j.out
#SBATCH --error=logs/dl_clean_lupus_%j.err
#
# Step 1/2 of the lupus pipeline: download from CELLxGENE Discover and
# clean. clean_lupus pulls clinical variables from a UCSF Box URL (HTTPS),
# so this step needs outbound internet — kept on cpu_p which has it.
# Preprocess is in scripts/preprocess_lupus.sbatch (submit with
# --dependency=afterok:<this jobid> for auto-chain).

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs data/lupus

CXG_ID="6e65a13e-2a45-4340-9442-bd2c41a01f17"
RAW=data/lupus/lupus.h5ad
CLEANED=data/lupus/lupus_cleaned.h5ad

if [ -f "$RAW" ]; then
    echo "===== raw exists, skipping download ($(ls -lh $RAW | awk '{print $5}')) ====="
else
    echo "===== download_from_cxg ($(date)) ====="
    ./target/executable/download_from_cxg/download_from_cxg \
        --input "$CXG_ID" \
        --output "$RAW"
    echo "downloaded: $(ls -lh $RAW | awk '{print $5}')"
fi

echo "===== clean_lupus ($(date)) ====="
./target/executable/clean_lupus/clean_lupus \
    --input "$RAW" \
    --output "$CLEANED" \
    --output_compression gzip

echo "===== done ($(date)) — cleaned: $(ls -lh $CLEANED | awk '{print $5}') ====="
