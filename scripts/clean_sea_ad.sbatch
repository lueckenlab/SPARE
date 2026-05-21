#!/bin/bash
#SBATCH --job-name=clean_seaad
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=600G
#SBATCH --time=12:00:00
#SBATCH --nice=10000
#SBATCH --output=logs/clean_sea_ad_%x_%j.out
#SBATCH --error=logs/clean_sea_ad_%x_%j.err
#
# Clean a SEA-AD Whole Taxonomy region via the executable runner.
# Usage: sbatch -J clean_sea_ad_mtg scripts/clean_sea_ad.sbatch mtg
#        sbatch -J clean_sea_ad_dlpfc scripts/clean_sea_ad.sbatch dlpfc
#
# ~1.4M nuclei x 35,483 genes; the CELLxGENE h5ad is ~50GB on disk and
# holds both raw (.raw) and normalized (.X) matrices (~56GB each in
# memory). The cleaner drops the normalized copy after the read, but peak
# (read + raw.to_adata + QC CSR rebuild) lands ~200GB, so 600G on a
# cpu_p big-memory node gives comfortable headroom.

set -euo pipefail

REGION="${1:?usage: clean_sea_ad.sbatch <mtg|dlpfc>}"
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
./target/executable/clean_sea_ad/clean_sea_ad \
    --input  "data/${DS}/${DS}.h5ad" \
    --output "data/${DS}/${DS}_cleaned.h5ad" \
    --output_compression gzip
echo "Wallclock: ${SECONDS}s"
echo "MaxRSS appears in sacct after the job ends."
