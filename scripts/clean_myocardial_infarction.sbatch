#!/bin/bash
#SBATCH --job-name=clean_mi
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --nice=10000
#SBATCH --output=logs/clean_myocardial_infarction_%j.out
#SBATCH --error=logs/clean_myocardial_infarction_%j.err
#
# Clean the Kuppe et al. 2022 myocardial infarction snRNA-seq atlas via
# the executable runner. 191,795 cells x 28,380 genes (~1.8GB on disk).
# The interactive node sandbox SIGKILLs the load at ~2GB RSS, so run it
# on a cpu_p node. 128GB is generous for a dataset this size (peak is a
# few GB: full .X + .raw + the QC-filter CSR rebuild).

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs

SECONDS=0
./target/executable/clean_myocardial_infarction/clean_myocardial_infarction \
    --input  data/myocardial_infarction/myocardial_infarction.h5ad \
    --output data/myocardial_infarction/myocardial_infarction_cleaned.h5ad \
    --output_compression gzip
echo "Wallclock: ${SECONDS}s"
echo "MaxRSS will appear in sacct after the job ends."
