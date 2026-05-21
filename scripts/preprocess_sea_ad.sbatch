#!/bin/bash
#SBATCH --job-name=pp_seaad
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=400G
#SBATCH --time=23:59:00
#SBATCH --nice=10000
#SBATCH --output=logs/pp_sea_ad_%x_%j.out
#SBATCH --error=logs/pp_sea_ad_%x_%j.err
#
# Preprocess a SEA-AD Whole Taxonomy region (~1.4M nuclei x 35,483 genes).
# Usage: sbatch -J pp_sea_ad_mtg   scripts/preprocess_sea_ad.sbatch mtg
#        sbatch -J pp_sea_ad_dlpfc scripts/preprocess_sea_ad.sbatch dlpfc
#
# HVG runs on the full gene set (~56GB raw load) before subsetting to
# 3000, then scVI/scANVI per batch covariate (donor_id, assay) + scPoli.
# imm_of_aging (3.76M cells) peaked ~325GB; at ~1.4M cells SEA-AD is
# smaller, so 400G + one GPU on a big node is ample. gpu_normal (24h)
# fits: scVI auto-scales epochs to cell count and scANVI is fixed at 20.
#
# Embeddings produced: X_pca, X_scVI_donor_id, X_scVI_assay,
# X_scANVI_donor_id, X_scANVI_assay, X_scpoli.

set -euo pipefail

REGION="${1:?usage: preprocess_sea_ad.sbatch <mtg|dlpfc>}"
case "$REGION" in
  mtg|dlpfc) ;;
  *) echo "REGION must be mtg or dlpfc" >&2; exit 2 ;;
esac
DS="sea_ad_${REGION}"

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

source /home/icb/vladimir.shitov/software/miniconda3/etc/profile.d/conda.sh
conda activate sample_representation_benchmark

mkdir -p logs

SECONDS=0
./target/executable/preprocess/preprocess \
    --input  "data/${DS}/${DS}_cleaned.h5ad" \
    --output "data/${DS}/${DS}_processed.h5ad" \
    --output_metadata "data/${DS}/${DS}_metadata.csv" \
    --sample_key "donor_id" \
    --cell_type_key "Subclass" \
    --batch_key "assay" \
    --batch_covariates "donor_id,assay" \
    --samples_metadata_cols "disease,sex,Neurotypical reference,Cognitive status,ADNC,Braak stage,Thal phase,CERAD score,APOE4 status,Age at death,Years of education,PMI,Lewy body disease pathology,LATE-NC stage,Microinfarct pathology,self_reported_ethnicity,development_stage,age,adnc_score,braak_score,cerad_score,thal_score" \
    --sample_size_threshold 100 \
    --output_compression gzip
echo "Wallclock: ${SECONDS}s"
