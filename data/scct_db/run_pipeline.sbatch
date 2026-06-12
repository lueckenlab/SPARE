#!/bin/bash
#SBATCH --job-name=scct_pipe
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/scct_pipe_%j.out
#SBATCH --error=logs/scct_pipe_%j.err
#
# Nextflow head process for the sample_representation pipeline over BOTH
# CD45 compartments (datasets_scct_db.yaml lists both). Child represent /
# aggregate / evaluate jobs are submitted to slurm by the `standard` profile.

set -euo pipefail
cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark
source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark
mkdir -p logs

RESUME_FLAG=""
if [ -d work_scct_db ]; then RESUME_FLAG="-resume"; fi

nextflow run target/nextflow/sample_representation/main.nf \
    -profile standard \
    -work-dir work_scct_db/ \
    --param_list src/workflows/sample_representation/datasets_scct_db.yaml \
    --method_params src/workflows/sample_representation/experiments_scct_db.yaml \
    --publish_dir output_scct_db/ \
    $RESUME_FLAG
echo "PIPELINE DONE"
