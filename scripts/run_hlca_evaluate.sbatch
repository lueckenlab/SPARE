#!/bin/bash
#SBATCH --job-name=eval_hlca
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/hlca_evaluate_%j.out
#SBATCH --error=logs/hlca_evaluate_%j.err

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs output_hlca_evaluate

./target/executable/evaluate/evaluate \
    --input data/hlca/hlca_representations.h5ad \
    --output_dir output_hlca_evaluate/ \
    --benchmark_schema data/hlca/benchmark_schema.json \
    --root_sample homosapiens_None_2023_None_sikkemalisa_002_d10_1101_2022_03_10_483747D322 \
    --trajectory_variable anatomical_region_ccf_score \
    --inverse_trajectory false \
    --figure_format png
