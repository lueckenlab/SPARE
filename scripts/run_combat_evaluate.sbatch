#!/bin/bash
#SBATCH --job-name=eval_combat
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/combat_evaluate_%j.out
#SBATCH --error=logs/combat_evaluate_%j.err

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs output_combat_evaluate

./target/executable/evaluate/evaluate \
    --input data/combat/combat_representations.h5ad \
    --output_dir output_combat_evaluate/ \
    --benchmark_schema data/combat/benchmark_schema.json \
    --root_sample H00052-Ha001E-PBGa \
    --trajectory_variable Outcome \
    --inverse_trajectory true \
    --figure_format png
