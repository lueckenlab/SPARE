#!/bin/bash
#SBATCH --job-name=spare_pilotgmvae
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --gres=gpu:1
#SBATCH --time=07:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/newmethods_combat_pilotgmvae_%j.out
#SBATCH --error=logs/newmethods_combat_pilotgmvae_%j.err
#
# Verification run of the pilot_gm_vae represent method on the FULL combat
# dataset (783k cells, 138 samples). Runs in the dedicated pilot_gm_vae env
# (see envs/README.md) and uses the GPU (cuda=1).

set -euo pipefail

REPO=/ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate pilot_gm_vae

IN=$REPO/data/combat/combat_processed.h5ad
OUT=$REPO/output_newmethods_combat
mkdir -p "$OUT"
# Run inside OUT so pilotgm's ./trained_models/ artifact does not pollute the repo.
cd "$OUT"

echo "===== pilot_gm_vae ($(date)) ====="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

"$REPO/target/executable/pilot_gm_vae/pilot_gm_vae" \
    --input "$IN" \
    --output "$OUT/pilot_gm_vae_pca.csv" \
    --sample_key scRNASeq_sample_ID \
    --cell_type_key Annotation_major_subset \
    --layer X_pca \
    --epochs 50 \
    --cuda 1

echo "===== done ($(date)) ====="
