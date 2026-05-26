#!/bin/bash
#SBATCH --job-name=eval_copd
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/copd_evaluate_%j.out
#SBATCH --error=logs/copd_evaluate_%j.err

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs output_copd_evaluate

./target/executable/evaluate/evaluate \
    --input data/copd/copd_representations.h5ad \
    --output_dir output_copd_evaluate/ \
    --benchmark_schema data/copd/benchmark_schema.json \
    --root_sample UA65_Batch1 \
    --trajectory_variable GOLD_Stage \
    --inverse_trajectory false \
    --figure_format png
