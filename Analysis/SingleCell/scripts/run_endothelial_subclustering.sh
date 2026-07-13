#!/bin/bash
# SBATCH orchestrator for endothelial (LSEC) + T-cell subclustering (WS1 / GLP1R).
# Chains: 313 (scVI, GPU) -> 314 (endothelial annotation, CPU)
#         313 (scVI, GPU) -> 315 (T-cell annotation, CPU)
#
# Usage: bash run_endothelial_subclustering.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INPUTS="${BASE}/Analysis/SingleCell/results_gpu_v2/mcp/inputs"
LOG_DIR="${SCRIPT_DIR}/logs/endothelial_subclustering"
mkdir -p "${LOG_DIR}"

ENDO_OUT="${BASE}/Analysis/SingleCell/results_gpu_v2/endothelial_subtypes"
TCELL_OUT="${BASE}/Analysis/SingleCell/results_gpu_v2/tcell_subtypes"
RESULTS="${BASE}/RNA-seq/results/glp1ra/lsec"
mkdir -p "${ENDO_OUT}" "${TCELL_OUT}" "${RESULTS}"

GPU_LD='export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/lib/python3.12/site-packages/nvidia/cu13/lib'

echo "=== Endothelial + T-cell subclustering (WS1) ==="
echo "Log dir: ${LOG_DIR}"

# ── 313 endothelial: scVI (GPU) ─────────────────────────────────────────
JOB_ENDO_SCVI=$(sbatch --parsable \
  --job-name=scvi --partition=gpu --gres=gpu:1 \
  --cpus-per-task=16 --mem=200G --time=48:00:00 \
  --output="${LOG_DIR}/313_endo_scvi_%j.out" \
  --error="${LOG_DIR}/313_endo_scvi_%j.err" \
  --wrap="${GPU_LD} && micromamba run -n rapids_singlecell python ${SCRIPT_DIR}/313_endothelial_subcluster_scvi.py \
      --input ${INPUTS}/atlas_cnmf_endothelial_cells.h5ad --outdir ${ENDO_OUT} --tag endothelial")
echo "313 endothelial scVI: ${JOB_ENDO_SCVI}"

# ── 314 endothelial: annotation + receptor query (CPU) ──────────────────
JOB_ENDO_ANNOT=$(sbatch --parsable \
  --dependency=afterok:${JOB_ENDO_SCVI} \
  --job-name=scanpy --partition=cpu \
  --cpus-per-task=8 --mem=128G --time=48:00:00 \
  --output="${LOG_DIR}/314_endo_annot_%j.out" \
  --error="${LOG_DIR}/314_endo_annot_%j.err" \
  --wrap="micromamba run -n spatial python ${SCRIPT_DIR}/314_endothelial_annotation.py \
      --outdir ${ENDO_OUT} --tag endothelial --results-dir ${RESULTS}")
echo "314 endothelial annotation: ${JOB_ENDO_ANNOT} (after ${JOB_ENDO_SCVI})"

# ── 313 T-cell: scVI (GPU) ──────────────────────────────────────────────
JOB_TCELL_SCVI=$(sbatch --parsable \
  --job-name=scvi --partition=gpu --gres=gpu:1 \
  --cpus-per-task=16 --mem=200G --time=48:00:00 \
  --output="${LOG_DIR}/313_tcell_scvi_%j.out" \
  --error="${LOG_DIR}/313_tcell_scvi_%j.err" \
  --wrap="${GPU_LD} && micromamba run -n rapids_singlecell python ${SCRIPT_DIR}/313_endothelial_subcluster_scvi.py \
      --input ${INPUTS}/atlas_cnmf_tcells.h5ad --outdir ${TCELL_OUT} --tag tcell")
echo "313 T-cell scVI: ${JOB_TCELL_SCVI}"

# ── 315 T-cell: annotation + GLP1R query (CPU) ──────────────────────────
JOB_TCELL_ANNOT=$(sbatch --parsable \
  --dependency=afterok:${JOB_TCELL_SCVI} \
  --job-name=scanpy --partition=cpu \
  --cpus-per-task=8 --mem=128G --time=48:00:00 \
  --output="${LOG_DIR}/315_tcell_annot_%j.out" \
  --error="${LOG_DIR}/315_tcell_annot_%j.err" \
  --wrap="micromamba run -n spatial python ${SCRIPT_DIR}/315_tcell_annotation.py \
      --outdir ${TCELL_OUT} --tag tcell --results-dir ${RESULTS}")
echo "315 T-cell annotation: ${JOB_TCELL_ANNOT} (after ${JOB_TCELL_SCVI})"

echo ""
echo "=== Submitted ==="
echo "313 endo scVI:   ${JOB_ENDO_SCVI}"
echo "314 endo annot:  ${JOB_ENDO_ANNOT}"
echo "313 tcell scVI:  ${JOB_TCELL_SCVI}"
echo "315 tcell annot: ${JOB_TCELL_ANNOT}"
