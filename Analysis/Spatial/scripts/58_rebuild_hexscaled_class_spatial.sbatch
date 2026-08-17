#!/usr/bin/env bash
#SBATCH --job-name=squidpy
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/hmsma_hexscaled_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/hmsma_hexscaled_%j.out
set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Rebuild the HMSMA AnnData objects with physically scaled hex coordinates, then
# recompute Moran's I and the matched-gene-set inference on the corrected graph.
#
# The existing HRA007511_starsolo build is NOT touched: the sealed program-arm
# candidate and the mask-sensitivity diagnostic read it directly, and the handoff
# forbids mutating a validated candidate. Those arms still carry values computed
# on the uncorrected graph and need their own versioned rebuild.
BUILD=Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled
RESULTS=Analysis/Spatial/results/hmsma_class_spatial_v3

micromamba run -n spatial python \
  Analysis/Spatial/scripts/53_build_hmsma_anndata.py \
  --out-dir "${BUILD}" || { echo "BUILD FAILED"; exit 1; }

# No --from-per-gene here: the cached Moran values were computed on the
# uncorrected graph and cannot be reused.
micromamba run -n spatial python \
  Analysis/Spatial/scripts/55_hmsma_class_spatial_autocorr.py \
  --in-dir "${BUILD}" \
  --out-dir "${RESULTS}" \
  --n-null 9999 \
  --seed 42 || { echo "ANALYSIS FAILED"; exit 1; }

echo "EXIT=0"
