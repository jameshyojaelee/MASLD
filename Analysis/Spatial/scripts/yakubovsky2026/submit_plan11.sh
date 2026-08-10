#!/usr/bin/env bash
set -euo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPTS="${BASE}/Analysis/Spatial/scripts/yakubovsky2026"
OUT="${BASE}/Analysis/Spatial/candidates/program-context-v2-candidate-2026-08-07/yakubovsky2026"
mkdir -p "${OUT}/logs"

adapter_job="$(sbatch --parsable \
  --output="${OUT}/logs/h5py_%j.out" \
  --error="${OUT}/logs/h5py_%j.err" \
  "${SCRIPTS}/run_adapter.sbatch")"
model_job="$(sbatch --parsable \
  --dependency="afterok:${adapter_job}" \
  --output="${OUT}/logs/bootstrap_%j.out" \
  --error="${OUT}/logs/bootstrap_%j.err" \
  "${SCRIPTS}/run_models.sbatch")"
validation_job="$(sbatch --parsable \
  --dependency="afterok:${model_job}" \
  --output="${OUT}/logs/python_%j.out" \
  --error="${OUT}/logs/python_%j.err" \
  "${SCRIPTS}/run_validate.sbatch")"

printf 'adapter_job=%s\nmodel_job=%s\nvalidation_job=%s\n' \
  "${adapter_job}" "${model_job}" "${validation_job}"

