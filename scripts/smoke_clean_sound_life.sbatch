#!/bin/bash
#SBATCH --job-name=clean_sl
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=900G
#SBATCH --time=8:00:00
#SBATCH --nice=10000
#SBATCH --output=logs/clean_sound_life_%j.out
#SBATCH --error=logs/clean_sound_life_%j.err
#
# Smoke run of src/clean_data/sound_life via the executable runner.
# sound_life is ~2.5x the nnz of imm_of_aging (22.9B vs 9.2B). The
# imm_of_aging smoke peaked at ~296GB (≈32 bytes/nnz, dominated by the
# QC filter's CSR rebuild), so sound_life is expected near ~730GB.
# 900GB lands on supercpu01/02 (2TB nodes in cpu_p). Compression off to
# avoid the gzip write buffer (re-enable once peak is confirmed).

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs

SECONDS=0
./target/executable/clean_sound_life/clean_sound_life \
    --input  data/sound_life/sound_life.h5ad \
    --output data/sound_life/sound_life_cleaned.h5ad
echo "Wallclock: ${SECONDS}s"
echo "MaxRSS will appear in sacct after the job ends."
