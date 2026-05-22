#!/bin/bash
#SBATCH --job-name=smoke_pilot_nf
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/smoke_pilot_nf_%j.out
#SBATCH --error=logs/smoke_pilot_nf_%j.err
#
# Smoke-test the pilot_gm_vae Nextflow env-activation path end-to-end on
# combat before the full matrix. This is the one untested piece: the
# component's beforeScript directive must source conda, activate the isolated
# pilot_gm_vae env, and reset LD_LIBRARY_PATH on a GPU node. Earlier pilot
# tests hit target/executable directly and never exercised this path.
#
# Nextflow itself is light (this job); the pilot_gm_vae job is submitted by
# Nextflow to gpu_p under the standard profile.

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs output_pilot_smoke

RESUME_FLAG=""
if [ -d work_pilot_smoke ] && [ -f .nextflow.log ]; then
    RESUME_FLAG="-resume"
fi

nextflow run target/nextflow/sample_representation/main.nf \
    -profile standard \
    -work-dir work_pilot_smoke/ \
    --param_list src/workflows/sample_representation/datasets_combat.yaml \
    --method_params src/workflows/sample_representation/experiments_pilot_smoke.yaml \
    --publish_dir output_pilot_smoke/ \
    $RESUME_FLAG

echo "===== nextflow exited; published outputs: ====="
ls -lhR output_pilot_smoke/ || true
