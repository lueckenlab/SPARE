#!/bin/bash
#SBATCH --job-name=spare_myocardial_infarction
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/myocardial_infarction_full_%j.out
#SBATCH --error=logs/myocardial_infarction_full_%j.err

set -euo pipefail

PROJECT_ROOT=/ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark
cd "$PROJECT_ROOT"

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs output_myocardial_infarction nf_runs/myocardial_infarction

RESUME_FLAG=""
if [ -d work_myocardial_infarction ]; then
    RESUME_FLAG="-resume"
fi

# Run nextflow from an isolated per-driver subdir so its `.nextflow/`
# (session locks, history) does not collide with the other concurrent
# represent drivers. All absolute paths back to the project resources.
cd nf_runs/myocardial_infarction
nextflow run "$PROJECT_ROOT/target/nextflow/sample_representation/main.nf" \
    -name "spare_myocardial_infarction_${SLURM_JOB_ID:-$$}" \
    -profile standard \
    -c "$PROJECT_ROOT/nextflow.config" \
    -work-dir "$PROJECT_ROOT/work_myocardial_infarction/" \
    --param_list "$PROJECT_ROOT/src/workflows/sample_representation/datasets_myocardial_infarction.yaml" \
    --method_params "$PROJECT_ROOT/src/workflows/sample_representation/experiments.yaml" \
    --publish_dir "$PROJECT_ROOT/output_myocardial_infarction/" \
    $RESUME_FLAG
