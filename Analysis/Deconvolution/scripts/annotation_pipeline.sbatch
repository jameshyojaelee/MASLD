#!/usr/bin/env bash
#SBATCH --job-name=annotation_pipeline
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --partition=cpu
#SBATCH --output=RNA-seq/deconvolution/logs/annotation_%j.out
#SBATCH --error=RNA-seq/deconvolution/logs/annotation_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13

# Use Conda environment's libraries instead of system
export LD_LIBRARY_PATH="RNA-seq/.mamba/music_deconv/lib:${LD_LIBRARY_PATH:-}"
# Prevent user local packages from interfering
export PYTHONNOUSERSITE=1

SPECIES=${SPECIES:-mouse}
MANIFEST=${MANIFEST:-Liver_Atlas/run_manifest.tsv}
MODEL=${MODEL:-Immune_All_Low}
OUTDIR=${OUTDIR:-RNA-seq/deconvolution/reference}

MICROMAMBA=RNA-seq/in-house_MCD_RNAseq/micromamba
ENV_PREFIX=RNA-seq/.mamba/music_deconv

mkdir -p "${OUTDIR}" RNA-seq/deconvolution/logs

REF_H5AD="${OUTDIR}/reference_${SPECIES}.h5ad"
CT_H5AD="${OUTDIR}/reference_${SPECIES}.celltypist.h5ad"
SV_H5AD="${OUTDIR}/reference_${SPECIES}.scanvi.h5ad"
FINAL_H5AD="${OUTDIR}/reference_${SPECIES}.consensus.h5ad"

echo "Starting annotation pipeline for ${SPECIES}"
echo "Manifest: ${MANIFEST}"
echo "Model: ${MODEL}"
echo "Output: ${OUTDIR}"

# Step 1: Build reference
echo "[Step 1/4] Building reference..."
"${MICROMAMBA}" run -p "${ENV_PREFIX}" python3 \
  RNA-seq/deconvolution/scripts/01_build_reference.py \
  --manifest "${MANIFEST}" \
  --species "${SPECIES}" \
  --output "${REF_H5AD}"

# Step 2: CellTypist annotation
echo "[Step 2/4] Running CellTypist..."
"${MICROMAMBA}" run -p "${ENV_PREFIX}" python3 \
  RNA-seq/deconvolution/scripts/02_celltypist.py \
  --input "${REF_H5AD}" \
  --output "${CT_H5AD}" \
  --model "${MODEL}" \
  --majority-voting

# Step 3: scANVI refinement
echo "[Step 3/4] Running scANVI..."
"${MICROMAMBA}" run -p "${ENV_PREFIX}" python3 \
  RNA-seq/deconvolution/scripts/03_scanvi.py \
  --input "${CT_H5AD}" \
  --output "${SV_H5AD}"

# Step 4: Consensus labels
echo "[Step 4/4] Computing consensus labels..."
"${MICROMAMBA}" run -p "${ENV_PREFIX}" python3 \
  RNA-seq/deconvolution/scripts/04_consensus_labels.py \
  --input "${SV_H5AD}" \
  --output "${FINAL_H5AD}"

echo "Annotation pipeline complete for ${SPECIES}!"
echo "Output: ${FINAL_H5AD}"
