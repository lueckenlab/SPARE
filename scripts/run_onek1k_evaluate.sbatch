#!/bin/bash
#SBATCH --job-name=eval_onek1k
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/onek1k_evaluate_%j.out
#SBATCH --error=logs/onek1k_evaluate_%j.err

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs output_onek1k_evaluate

./target/executable/evaluate/evaluate \
    --input data/onek1k/onek1k_representations.h5ad \
    --output_dir output_onek1k_evaluate/ \
    --benchmark_schema data/onek1k/benchmark_schema.json \
    --root_sample 177_178 \
    --trajectory_variable age \
    --inverse_trajectory false \
    --figure_format png
