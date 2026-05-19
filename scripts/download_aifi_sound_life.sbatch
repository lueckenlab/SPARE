#!/bin/bash
#SBATCH --job-name=aifi_dl_sl
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --array=0-7%4
#SBATCH --output=logs/aifi_dl_sl_%A_%a.log
#SBATCH --error=logs/aifi_dl_sl_%A_%a.log
#SBATCH --nice=10000
#
# Downloads the 8 per-cell-type h5ad files of the AIFI Sound Life
# vaccine cohort into data/sound_life/. Adapted from
# patpy-vaccine-prediction/scripts/download_aifi_sound_life.sbatch.
#
# Reproducibility-only: in normal SPARE runs the data is already on
# disk and symlinked from data/sound_life/. After running this you
# still need to concatenate the per-cell-type files into a single
# sound_life.h5ad (see patpy-vaccine-prediction/scripts/concat_aifi.py).

set -euo pipefail

cd /ictstr01/groups/luckylab/workspace/vladimir.shitov/2023_05_patient_representation_benchmark/reproducibility/pat_rep_benchmark

DATA_DIR="data/sound_life"
MANIFEST="${DATA_DIR}/urls.txt"

mkdir -p logs "${DATA_DIR}"

LINE_NO=$(( SLURM_ARRAY_TASK_ID + 1 ))
LINE=$(sed -n "${LINE_NO}p" "${MANIFEST}")
FNAME=$(awk -F'\t' '{print $1}' <<<"${LINE}")
URL=$(awk   -F'\t' '{print $2}' <<<"${LINE}")
EXPECTED=$(awk -F'\t' '{print $3}' <<<"${LINE}")
DEST="${DATA_DIR}/${FNAME}"

echo "=== $(date) === [task ${SLURM_ARRAY_TASK_ID}] ${FNAME}"
echo "URL:      ${URL}"
echo "Dest:     ${DEST}"
echo "Expected: ${EXPECTED} bytes"

MIN_OK=$(( EXPECTED * 80 / 100 ))
if [[ -f "${DEST}" ]]; then
    HAVE=$(stat -c%s "${DEST}")
    if (( HAVE >= MIN_OK )); then
        echo "Already complete (${HAVE} bytes, >=${MIN_OK}). Skipping."
        exit 0
    fi
    echo "Partial file present (${HAVE} bytes < ${MIN_OK}); restarting (no range support)."
fi

PARTIAL="${DEST}.partial"
rm -f "${PARTIAL}"

curl --fail --location --silent --show-error \
     --http1.1 \
     --speed-time 60 --speed-limit 1024 \
     --retry 10 --retry-delay 30 --retry-all-errors \
     -o "${PARTIAL}" "${URL}"

HAVE=$(stat -c%s "${PARTIAL}")
echo "Got ${HAVE} bytes (expected ~${EXPECTED}, min ${MIN_OK})."
if (( HAVE < MIN_OK )); then
    echo "ERROR: file truncated. Leaving .partial for inspection."
    exit 1
fi

~/software/miniconda3/envs/sample_representation_benchmark/bin/python -c "import h5py,sys;sys.exit(0 if h5py.is_hdf5('${PARTIAL}') else 1)" || {
    echo "ERROR: not a valid HDF5 file."
    exit 1
}

mv -f "${PARTIAL}" "${DEST}"
echo "=== $(date) === Done: ${DEST}"
