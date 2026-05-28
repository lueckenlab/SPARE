#!/bin/bash
#SBATCH --job-name=spare_onek1k
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/onek1k_full_%j.out
#SBATCH --error=logs/onek1k_full_%j.err

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs output_onek1k

RESUME_FLAG=""
if [ -d work_onek1k ] && [ -f .nextflow.log ]; then
    RESUME_FLAG="-resume"
fi

nextflow run target/nextflow/sample_representation/main.nf \
    -profile standard \
    -work-dir work_onek1k/ \
    --param_list src/workflows/sample_representation/datasets_onek1k.yaml \
    --method_params src/workflows/sample_representation/experiments.yaml \
    --publish_dir output_onek1k/ \
    $RESUME_FLAG
