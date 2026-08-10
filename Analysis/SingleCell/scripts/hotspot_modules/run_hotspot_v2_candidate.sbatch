#!/bin/bash
#SBATCH --job-name=hotspot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/logs/hotspot_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/logs/hotspot_%j.err

set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
RELEASE_ID="program-context-v2-candidate-2026-08-07"
OUT_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/${RELEASE_ID}/hotspot"
RNASEQ_RSCRIPT="${RNASEQ_RSCRIPT:-/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript}"

# The log directory must exist before sbatch opens --output/--error. Prepare it
# before submission with:
#   mkdir -p Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/logs
test -d "${OUT_ROOT}/logs" || {
  echo "Missing pre-created SLURM log directory: ${OUT_ROOT}/logs" >&2
  exit 2
}

cd "${PROJECT_ROOT}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export HOTSPOT_V2_RELEASE_ID="${RELEASE_ID}"
test -x "${RNASEQ_RSCRIPT}" || {
  echo "Missing executable rnaseq Rscript: ${RNASEQ_RSCRIPT}" >&2
  exit 2
}

case "${HOTSPOT_V2_CHECK_ONLY:-false}" in
  true|TRUE|t|T|1|yes|YES|y|Y)
    "${RNASEQ_RSCRIPT}" \
      Analysis/SingleCell/scripts/hotspot_modules/512_hotspot_v2_donor_refit.R \
      --check-only
    "${RNASEQ_RSCRIPT}" \
      Analysis/SingleCell/scripts/hotspot_modules/513_freeze_hotspot_v2_registry.R \
      --check-only
    echo "Hotspot v2 check-only preflight completed; no candidate outputs written"
    exit 0
    ;;
esac

"${RNASEQ_RSCRIPT}" \
  Analysis/SingleCell/scripts/hotspot_modules/512_hotspot_v2_donor_refit.R

"${RNASEQ_RSCRIPT}" \
  Analysis/SingleCell/scripts/hotspot_modules/513_freeze_hotspot_v2_registry.R

"${RNASEQ_RSCRIPT}" \
  Analysis/SingleCell/scripts/hotspot_modules/514_validate_hotspot_v2_registry.R

test -s "${OUT_ROOT}/validation_status.tsv"
test -s "${OUT_ROOT}/READY"
echo "Hotspot v2 candidate frozen and independently validated under ${OUT_ROOT}"
