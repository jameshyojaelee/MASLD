#!/usr/bin/env bash
#SBATCH --job-name=masld_cl_scope
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/scope-2026-08-13T231950Z/slurm.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/scope-2026-08-13T231950Z/slurm.log

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RUN_ROOT="${PROJECT_ROOT}/Analysis/SingleCell/continual_integration/candidates/scope-2026-08-13T231950Z"
cd "${PROJECT_ROOT}"
trap 'status=$?; if [[ ${status} -eq 0 ]]; then touch "${RUN_ROOT}/SCOPE_COMPLETE"; else printf "%s\n" "${status}" > "${RUN_ROOT}/FAILED.exit_code"; fi' EXIT
sha256sum docs/STATUS.md docs/ROADMAP.md scripts/manuscript/validate_resource_scope.py \
  > "${RUN_ROOT}/input_sha256.txt"
python3 scripts/manuscript/validate_resource_scope.py \
  > "${RUN_ROOT}/validation.txt"
