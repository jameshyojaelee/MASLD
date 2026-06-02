#!/usr/bin/env bash
#SBATCH --job-name=deconv_GSE213621
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/GSE213621_deconv_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/GSE213621_deconv_%j.err

# GSE213621 Deconvolution: MuSiC + BayesPrism (InstaPrism)
# 367 samples — largest human MASLD cohort
# Previously missing from deconvolution pipeline (discovered during Phase A monitoring)
#
# Steps:
#   1. Convert SRA metadata to TSV
#   2. Prepare bulk counts (Ensembl → symbol mapping)
#   3. Run MuSiC deconvolution
#   4. Run QC on MuSiC results
#   5. Run BayesPrism (InstaPrism) deconvolution
#
# NOTE: Uses R directly (not Rscript) because the Rscript binary in music_deconv
# env has a hard-coded R_HOME from the original install location. R itself works fine.
# This matches the pattern used in run_bayesprism_missing_human.sbatch.

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV_PREFIX="${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv"
SCRIPTS_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/scripts"

REF_HUMAN="${PROJECT_ROOT}/Analysis/Deconvolution/reference/reference_human_sce.rds"
BULK_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/bulk"
RESULTS_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/results"
META_BASE="${PROJECT_ROOT}/Analysis/Deconvolution/metadata"

DATASET="GSE213621"
SPECIES="human"
COUNTS_PATH="${PROJECT_ROOT}/RNA-seq/Human/Patient_Cohorts/results/${DATASET}/counts/featurecounts/gene_counts.txt"
SRA_METADATA="${PROJECT_ROOT}/RNA-seq/Human/Patient_Cohorts/pipelines/custom/${DATASET}/metadata/SraRunTable.csv"

# Use R binary directly (not Rscript — known env issue)
R_BIN="${ENV_PREFIX}/bin/R"
export R_LIBS_USER=""
export R_LIBS_SITE="${ENV_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

echo "=============================="
echo "Deconvolution: ${DATASET}"
echo "Host: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-16}"
echo "Started: $(date)"
echo "R binary: ${R_BIN}"
echo "=============================="

# Verify inputs
for f in "${COUNTS_PATH}" "${SRA_METADATA}" "${REF_HUMAN}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: Missing input file: ${f}" >&2
    exit 1
  fi
done

echo "All input files verified."

# Paths for intermediate and output files
METADATA_TSV="${META_BASE}/${DATASET}_sra_metadata.tsv"
BULK_DIR="${BULK_BASE}/${DATASET}"
RESULT_DIR="${RESULTS_BASE}/${DATASET}"
QC_DIR="${RESULT_DIR}/qc"
COUNTS_TSV="${BULK_DIR}/${DATASET}_counts.tsv"
BULK_META_TSV="${BULK_DIR}/${DATASET}_metadata.tsv"

mkdir -p "${BULK_DIR}" "${RESULT_DIR}" "${QC_DIR}"

# ── Step 1: Convert SRA metadata ──
echo ""
echo "=== Step 1: Convert SRA metadata ==="
if [[ -f "${METADATA_TSV}" ]]; then
  echo "  Metadata TSV already exists, skipping conversion."
else
  python3 "${SCRIPTS_DIR}/00_convert_sra_metadata.py" \
    --input "${SRA_METADATA}" \
    --output "${METADATA_TSV}"
  echo "  Done. Output: ${METADATA_TSV}"
fi

# Verify metadata was created
if [[ ! -f "${METADATA_TSV}" ]]; then
  echo "ERROR: Metadata conversion failed." >&2
  exit 1
fi
echo "  Metadata samples: $(( $(wc -l < "${METADATA_TSV}") - 1 ))"

# ── Step 2: Prepare bulk counts ──
echo ""
echo "=== Step 2: Prepare bulk counts ==="
if [[ -f "${COUNTS_TSV}" && -f "${BULK_META_TSV}" ]]; then
  echo "  Bulk counts already prepared, skipping."
else
  "${R_BIN}" --no-save --no-restore --slave \
    --file="${SCRIPTS_DIR}/06_prepare_bulk_counts.R" \
    --args \
    "${COUNTS_PATH}" \
    "${METADATA_TSV}" \
    "${SPECIES}" \
    "${BULK_DIR}" \
    "${DATASET}"
  echo "  Done."
fi

# Verify counts prepared
if [[ ! -f "${COUNTS_TSV}" ]]; then
  echo "ERROR: Bulk count preparation failed — ${COUNTS_TSV} not found." >&2
  exit 1
fi
echo "  Counts dimensions: $(head -1 "${COUNTS_TSV}" | tr '\t' '\n' | wc -l) columns"

# ── Step 3: Run MuSiC ──
echo ""
echo "=== Step 3: Run MuSiC deconvolution ==="
MUSIC_OUT="${RESULT_DIR}/${DATASET}_music_prop_weighted.tsv"
if [[ -f "${MUSIC_OUT}" ]]; then
  echo "  MuSiC output already exists, skipping."
else
  "${R_BIN}" --no-save --no-restore --slave \
    --file="${SCRIPTS_DIR}/07_run_music.R" \
    --args \
    "${REF_HUMAN}" \
    "${COUNTS_TSV}" \
    "${BULK_META_TSV}" \
    "${RESULT_DIR}" \
    "${DATASET}"
  echo "  Done."
fi

# Verify MuSiC output
if [[ ! -f "${MUSIC_OUT}" ]]; then
  echo "ERROR: MuSiC deconvolution failed — ${MUSIC_OUT} not found." >&2
  exit 1
fi
echo "  MuSiC output samples: $(( $(wc -l < "${MUSIC_OUT}") - 1 ))"

# ── Step 4: QC on MuSiC results ──
echo ""
echo "=== Step 4: QC on MuSiC results ==="
"${R_BIN}" --no-save --no-restore --slave \
  --file="${SCRIPTS_DIR}/08_qc_music.R" \
  --args \
  "${MUSIC_OUT}" \
  "${BULK_META_TSV}" \
  "${QC_DIR}" \
  "${DATASET}"
echo "  Done."

# ── Step 5: Run BayesPrism (InstaPrism) ──
echo ""
echo "=== Step 5: Run BayesPrism (InstaPrism) deconvolution ==="
BP_OUT="${RESULT_DIR}/${DATASET}_bayesprism_proportions.tsv"
if [[ -f "${BP_OUT}" ]]; then
  echo "  BayesPrism output already exists, skipping."
else
  "${R_BIN}" --no-save --no-restore --slave \
    --file="${SCRIPTS_DIR}/12_run_bayesprism.R" \
    --args \
    "${REF_HUMAN}" \
    "${COUNTS_TSV}" \
    "${RESULT_DIR}" \
    "${DATASET}" \
    "celltype"
  echo "  Done."
fi

# Verify BayesPrism output
if [[ ! -f "${BP_OUT}" ]]; then
  echo "ERROR: BayesPrism deconvolution failed — ${BP_OUT} not found." >&2
  exit 1
fi
echo "  BayesPrism output samples: $(( $(wc -l < "${BP_OUT}") - 1 ))"

# ── Summary ──
echo ""
echo "=============================="
echo "All steps completed successfully for ${DATASET}."
echo "Finished: $(date)"
echo ""
echo "Outputs:"
echo "  Metadata: ${METADATA_TSV}"
echo "  Counts:   ${COUNTS_TSV}"
echo "  MuSiC:    ${MUSIC_OUT}"
echo "  BayesPrism: ${BP_OUT}"
echo "  QC:       ${QC_DIR}/"
echo ""
ls -lh "${RESULT_DIR}/"
echo "=============================="
