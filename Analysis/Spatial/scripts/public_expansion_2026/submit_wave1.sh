#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
SCRIPT_DIR="${PROJECT_ROOT}/Analysis/Spatial/scripts/public_expansion_2026"

YAK_LEAN_ID="$(sbatch --parsable "${SCRIPT_DIR}/run_yakubovsky_lean.sbatch")"
YAK_VISIUM_ID="$(sbatch --parsable "${SCRIPT_DIR}/run_yakubovsky_visium.sbatch")"
GEO_ID="$(sbatch --parsable "${SCRIPT_DIR}/run_gse287826.sbatch")"

printf 'task_id\tjob_id\n'
printf 'SP-AQ-YAK-LEAN\t%s\n' "${YAK_LEAN_ID}"
printf 'SP-AQ-YAK-VISIUM\t%s\n' "${YAK_VISIUM_ID}"
printf 'SP-AQ-GEO\t%s\n' "${GEO_ID}"
