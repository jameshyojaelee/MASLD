#!/bin/bash
# SBATCH orchestrator for hepatocyte subclustering pipeline
# Chains: 308 (GPU) → 309 (CPU) → 310a/b/c (CPU, parallel)
# 311 is NOT auto-chained — run manually after label review.
#
# Usage: bash run_hepatocyte_subclustering.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs/hepatocyte_subclustering"
mkdir -p "${LOG_DIR}"

echo "=== Hepatocyte Subclustering Pipeline ==="
echo "Log dir: ${LOG_DIR}"

# ── 308: scVI training (GPU) ────────────────────────────────────────────
JOB_308=$(sbatch --parsable \
  --job-name=hep_308_scvi \
  --partition=gpu \
  --gres=gpu:1 \
  --cpus-per-task=16 \
  --mem=200G \
  --time=48:00:00 \
  --output="${LOG_DIR}/308_scvi_%j.out" \
  --error="${LOG_DIR}/308_scvi_%j.err" \
  --wrap="export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/lib/python3.12/site-packages/nvidia/cu13/lib && micromamba run -n rapids_singlecell python ${SCRIPT_DIR}/308_hepatocyte_subcluster_scvi.py")
echo "308 submitted: ${JOB_308}"

# ── 309: Annotation (CPU, bigmem for 657K cells) ────────────────────────
JOB_309=$(sbatch --parsable \
  --dependency=afterok:${JOB_308} \
  --job-name=hep_309_annot \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=128G \
  --time=48:00:00 \
  --output="${LOG_DIR}/309_annot_%j.out" \
  --error="${LOG_DIR}/309_annot_%j.err" \
  --wrap="micromamba run -n spatial python ${SCRIPT_DIR}/309_hepatocyte_subcluster_annotation.py")
echo "309 submitted: ${JOB_309} (depends on ${JOB_308})"

# ── 310a: Bulk cross-modal (CPU) ────────────────────────────────────────
JOB_310A=$(sbatch --parsable \
  --dependency=afterok:${JOB_309} \
  --job-name=hep_310a_bulk \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=128G \
  --time=48:00:00 \
  --output="${LOG_DIR}/310a_bulk_%j.out" \
  --error="${LOG_DIR}/310a_bulk_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${SCRIPT_DIR}/310a_hepatocyte_crossmodal_bulk.R")
echo "310a submitted: ${JOB_310A} (depends on ${JOB_309})"

# ── 310b: Spatial + GWAS-ATAC (CPU) ─────────────────────────────────────
JOB_310B=$(sbatch --parsable \
  --dependency=afterok:${JOB_309} \
  --job-name=hep_310b_spatial \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=128G \
  --time=48:00:00 \
  --output="${LOG_DIR}/310b_spatial_%j.out" \
  --error="${LOG_DIR}/310b_spatial_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${SCRIPT_DIR}/310b_hepatocyte_crossmodal_spatial_gwas.R")
echo "310b submitted: ${JOB_310B} (depends on ${JOB_309})"

# ── 310c: Plasma (CPU) ──────────────────────────────────────────────────
JOB_310C=$(sbatch --parsable \
  --dependency=afterok:${JOB_309} \
  --job-name=hep_310c_plasma \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=64G \
  --time=48:00:00 \
  --output="${LOG_DIR}/310c_plasma_%j.out" \
  --error="${LOG_DIR}/310c_plasma_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${SCRIPT_DIR}/310c_hepatocyte_crossmodal_plasma.R")
echo "310c submitted: ${JOB_310C} (depends on ${JOB_309})"

echo ""
echo "=== Pipeline submitted ==="
echo "308 (scVI):   ${JOB_308}"
echo "309 (annot):  ${JOB_309}"
echo "310a (bulk):  ${JOB_310A}"
echo "310b (spat):  ${JOB_310B}"
echo "310c (plas):  ${JOB_310C}"
echo ""
echo "After 310a/b/c complete:"
echo "  1. Review hepatocyte_subtypes/subtype_disease_enrichment.csv"
echo "  2. Edit labels in hepatocyte_subtypes/hepatocyte_subtype_metadata.csv"
echo "  3. Run: sbatch ${SCRIPT_DIR}/run_hepatocyte_figures.sh"
