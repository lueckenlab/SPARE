#!/bin/bash
#SBATCH --job-name=spare_hlca
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/hlca_full_%j.out
#SBATCH --error=logs/hlca_full_%j.err

set -euo pipefail

PROJECT_ROOT=/ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark
cd "$PROJECT_ROOT"

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs output_hlca nf_runs/hlca

RESUME_FLAG=""
if [ -d work_hlca ]; then
    RESUME_FLAG="-resume"
fi

# Run nextflow from an isolated per-driver subdir so its `.nextflow/`
# (session locks, history) does not collide with the other concurrent
# represent drivers. All absolute paths back to the project resources.
cd nf_runs/hlca
nextflow run "$PROJECT_ROOT/target/nextflow/sample_representation/main.nf" \
    -name "spare_hlca_${SLURM_JOB_ID:-$$}" \
    -profile standard \
    -c "$PROJECT_ROOT/nextflow.config" \
    -work-dir "$PROJECT_ROOT/work_hlca/" \
    --param_list "$PROJECT_ROOT/src/workflows/sample_representation/datasets_hlca.yaml" \
    --method_params "$PROJECT_ROOT/src/workflows/sample_representation/experiments.yaml" \
    --publish_dir "$PROJECT_ROOT/output_hlca/" \
    $RESUME_FLAG
