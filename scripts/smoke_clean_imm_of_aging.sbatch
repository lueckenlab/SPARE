#!/bin/bash
#SBATCH --job-name=clean_ioa
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --time=6:00:00
#SBATCH --nice=10000
#SBATCH --output=logs/clean_imm_of_aging_%j.out
#SBATCH --error=logs/clean_imm_of_aging_%j.err
#
# Smoke run of src/clean_data/imm_of_aging via the executable runner.
# Drives the cleaning step end-to-end on the concat h5ad and writes
# data/imm_of_aging/imm_of_aging_cleaned.h5ad.
#
# Memory budget: peak observed at 249GB on the first attempt (OOM in
# step 36529206). The X load + int32 cast + layers["X_raw_counts"]
# copy + gzip write together exceed 250GB. Bumping to 400GB; resize
# after observing the new high-water mark.

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs

# Skip gzip on the smoke run — the previous attempt (36529684) OOM'd at
# 400GB during the gzip-buffered write. Uncompressed h5ad will be ~150GB
# on disk; we'll re-enable compression once we know the actual peak.
SECONDS=0
./target/executable/clean_imm_of_aging/clean_imm_of_aging \
    --input  data/imm_of_aging/imm_of_aging.h5ad \
    --output data/imm_of_aging/imm_of_aging_cleaned.h5ad
echo "Wallclock: ${SECONDS}s"
echo "MaxRSS will appear in sacct after the job ends."
