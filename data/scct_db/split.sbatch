#!/bin/bash
#SBATCH --job-name=scct_split
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --nice=10000
#SBATCH --output=logs/scct_split_%j.out
#SBATCH --error=logs/scct_split_%j.err

set -euo pipefail
cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark
source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark
python data/scct_db/split_by_sort.py
