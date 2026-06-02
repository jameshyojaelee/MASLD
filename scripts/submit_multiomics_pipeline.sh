#!/usr/bin/env bash
# Master orchestrator for multi-omics processing pipeline.
# Submits all 4 phases with SLURM dependency chains.
#
# Phase 1: CellRanger (GSE136103) — cpu, array 1-51
# Phase 2: scVI Integration        — bigmem, depends on Phase 1
# Phase 3: Proteomics Processing   — cpu, independent
# Phase 4: Proteo-Transcriptomic   — cpu, depends on Phase 2 + 3
#
# Usage: bash scripts/submit_multiomics_pipeline.sh

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJECT_ROOT}"

echo "============================================"
echo "  Multi-Omics Processing Pipeline Launcher"
echo "============================================"
echo ""

# Verify prerequisites
echo "Checking prerequisites..."
CELLRANGER="${HOME}/cellranger-10.0.0/cellranger"
if [[ ! -x "${CELLRANGER}" ]]; then
  echo "ERROR: CellRanger 10.0.0 not found at ${CELLRANGER}" >&2
  echo "  Run: bash scripts/setup/download_cellranger10.sh" >&2
  exit 1
fi
echo "  CellRanger 10.0.0: OK"

MANIFEST="data/GSE136103/run_manifest.tsv"
if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: Manifest not found: ${MANIFEST}" >&2
  echo "  Run: python3 data/GSE136103/scripts/build_manifest.py" >&2
  exit 1
fi
N_RUNS=$(awk -F'\t' 'NR>1 && $7!="1"' "${MANIFEST}" | wc -l)
echo "  GSE136103 manifest: ${N_RUNS} runs to process"

if [[ "${N_RUNS}" -eq 0 ]]; then
  echo "  WARNING: All CellRanger runs already complete. Skipping Phase 1."
  echo ""
fi

# Create log directories
mkdir -p data/GSE136103/logs
mkdir -p Analysis/SingleCell/logs
mkdir -p Analysis/Proteomics/logs

# Phase 1: CellRanger
echo ""
echo "--- Phase 1: CellRanger (GSE136103) ---"
if [[ "${N_RUNS}" -gt 0 ]]; then
  JOB_CR=$(sbatch --parsable --array=1-${N_RUNS}%8 \
    data/GSE136103/scripts/cellranger_array.sbatch)
  echo "  Submitted: job ${JOB_CR} (array 1-${N_RUNS}, 8 concurrent)"

  # Phase 1b: Scanpy QC (after cellranger)
  JOB_QC=$(sbatch --parsable --dependency=afterany:${JOB_CR} \
    --array=1-${N_RUNS}%16 \
    data/GSE136103/scripts/scanpy_qc.sbatch)
  echo "  Phase 1b (Scanpy QC): job ${JOB_QC}"

  PHASE1_DEP="--dependency=afterok:${JOB_QC}"
else
  PHASE1_DEP=""
fi

# Phase 2: scVI Integration
echo ""
echo "--- Phase 2: scVI Integration ---"
JOB_SCVI=$(sbatch --parsable ${PHASE1_DEP} \
  Analysis/SingleCell/scripts/run_scvi.sbatch)
echo "  Submitted: job ${JOB_SCVI}"

# Phase 3: Proteomics (independent)
echo ""
echo "--- Phase 3: Proteomics Processing ---"
JOB_PROT=$(sbatch --parsable \
  Analysis/Proteomics/scripts/run_proteomics.sbatch)
echo "  Submitted: job ${JOB_PROT}"

# Phase 4: Integration (after Phase 2 + 3)
echo ""
echo "--- Phase 4: Proteo-Transcriptomic Integration ---"
JOB_INT=$(sbatch --parsable --dependency=afterok:${JOB_SCVI}:${JOB_PROT} \
  Analysis/Proteomics/scripts/run_integration.sbatch)
echo "  Submitted: job ${JOB_INT}"

echo ""
echo "============================================"
echo "  All phases submitted!"
echo "============================================"
echo ""
echo "  Phase 1 (CellRanger):  ${JOB_CR:-SKIPPED}"
echo "  Phase 1b (Scanpy QC):  ${JOB_QC:-SKIPPED}"
echo "  Phase 2 (scVI):        ${JOB_SCVI}"
echo "  Phase 3 (Proteomics):  ${JOB_PROT}"
echo "  Phase 4 (Integration): ${JOB_INT}"
echo ""
echo "  Monitor: squeue -u ${USER}"
echo "  Logs:    data/GSE136103/logs/, Analysis/*/logs/"
