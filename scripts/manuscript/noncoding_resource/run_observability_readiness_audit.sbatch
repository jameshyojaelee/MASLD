#!/bin/bash
#SBATCH --job-name=ncobsaud
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncobsaud_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncobsaud_%j.err

set -euo pipefail
umask 027

: "${NONCODING_OBSERVABILITY_OUTPUT:?Set NONCODING_OBSERVABILITY_OUTPUT to a new output directory}"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PYTHON="/usr/bin/python3"
PRODUCER="${PROJECT_ROOT}/scripts/manuscript/noncoding_resource/audit_observability_readiness.py"
BIOTYPE_ROOT="${PROJECT_ROOT}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/biotype_association/run_gate19751171"
BIOTYPE_REVIEW="${PROJECT_ROOT}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/biotype_association/review_gate19751171_job19751607/VALIDATED.json"
EQTL_READY="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/genetics_context/observability_readiness.tsv"
SPATIAL_AXIS="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/coverage/gene_axis_audit.tsv"
SC_UNIVERSE="${PROJECT_ROOT}/Analysis/SingleCell/results_gpu_v2/mcp/inputs/dialogue_pseudobulk/gene_universe.tsv"

test ! -e "${NONCODING_OBSERVABILITY_OUTPUT}"
"${PYTHON}" "${PRODUCER}" \
  --biotype-counts "${BIOTYPE_ROOT}/bulk_deg_counts_by_biotype.tsv" \
  --biotype-review "${BIOTYPE_REVIEW}" \
  --eqtl-readiness "${EQTL_READY}" \
  --spatial-axis-audit "${SPATIAL_AXIS}" \
  --singlecell-gene-universe "${SC_UNIVERSE}" \
  --output "${NONCODING_OBSERVABILITY_OUTPUT}"

test -f "${NONCODING_OBSERVABILITY_OUTPUT}/VALIDATED.json"
grep -q '"scientific_gate": "skipped_incomplete_assay_matrix"' \
  "${NONCODING_OBSERVABILITY_OUTPUT}/VALIDATED.json"
