#!/bin/bash
#SBATCH --job-name=spare_stephenson
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/stephenson_full_%j.out
#SBATCH --error=logs/stephenson_full_%j.err

set -euo pipefail

PROJECT_ROOT=/ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark
cd "$PROJECT_ROOT"

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs output_stephenson nf_runs/stephenson

RESUME_FLAG=""
if [ -d work_stephenson ]; then
    RESUME_FLAG="-resume"
fi

# Run nextflow from an isolated per-driver subdir so its `.nextflow/`
# (session locks, history) does not collide with the other concurrent
# represent drivers. All absolute paths back to the project resources.
cd nf_runs/stephenson
nextflow run "$PROJECT_ROOT/target/nextflow/sample_representation/main.nf" \
    -name "spare_stephenson_${SLURM_JOB_ID:-$$}" \
    -profile standard \
    -c "$PROJECT_ROOT/nextflow.config" \
    -work-dir "$PROJECT_ROOT/work_stephenson/" \
    --param_list "$PROJECT_ROOT/src/workflows/sample_representation/datasets_stephenson.yaml" \
    --method_params "$PROJECT_ROOT/src/workflows/sample_representation/experiments.yaml" \
    --publish_dir "$PROJECT_ROOT/output_stephenson/" \
    $RESUME_FLAG
