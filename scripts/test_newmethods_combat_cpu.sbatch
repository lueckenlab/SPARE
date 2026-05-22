#!/bin/bash
#SBATCH --job-name=spare_newmethods_cpu
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=07:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/newmethods_combat_cpu_%j.out
#SBATCH --error=logs/newmethods_combat_cpu_%j.err
#
# Verification run of the two new shared-env represent methods (diffusionemd,
# phemd) on the FULL combat dataset (783k cells, 138 samples).

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

IN=data/combat/combat_processed.h5ad
OUT=output_newmethods_combat
mkdir -p "$OUT"

if [ -f "$OUT/diffusionemd_pca.csv" ]; then
    echo "===== diffusionemd: output exists, skipping ====="
else
    echo "===== diffusionemd ($(date)) ====="
    target/executable/diffusionemd/diffusionemd \
        --input "$IN" \
        --output "$OUT/diffusionemd_pca.csv" \
        --sample_key scRNASeq_sample_ID \
        --cell_type_key Annotation_major_subset \
        --layer X_pca
fi

echo "===== phemd ($(date)) ====="
target/executable/phemd/phemd \
    --input "$IN" \
    --output "$OUT/phemd_pca.csv" \
    --sample_key scRNASeq_sample_ID \
    --cell_type_key Annotation_major_subset \
    --layer X_pca \
    --subset_n_obs 2000 \
    --n_jobs 16

echo "===== done ($(date)) ====="
