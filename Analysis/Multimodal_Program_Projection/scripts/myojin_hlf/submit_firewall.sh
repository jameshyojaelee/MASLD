#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPT_DIR="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/myojin_hlf"

cd "${PROJECT_ROOT}"

custodian_job="$(sbatch --parsable "${SCRIPT_DIR}/run_custodian.sbatch")"
depmap_job="$(sbatch --parsable "${SCRIPT_DIR}/run_depmap.sbatch")"
freeze_job="$(sbatch --parsable --dependency="afterok:${custodian_job}:${depmap_job}" "${SCRIPT_DIR}/run_freeze.sbatch")"
printf 'custodian_job\t%s\ndepmap_job\t%s\nfreeze_job\t%s\n' "${custodian_job}" "${depmap_job}" "${freeze_job}"
