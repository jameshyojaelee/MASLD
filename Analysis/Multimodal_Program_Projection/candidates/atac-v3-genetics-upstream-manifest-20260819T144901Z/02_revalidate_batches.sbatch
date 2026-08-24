#!/bin/bash
#SBATCH --job-name=atacreval
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/atac-v3-genetics-upstream-manifest-20260819T144901Z/revalidate_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/atac-v3-genetics-upstream-manifest-20260819T144901Z/revalidate_%j.log

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
GEN=${ROOT}/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/genetics
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
eval "$(${MAMBA} shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

revalidate () {
  local batch=$1 pending=$2
  local final=${GEN}/replay_execution/batch_${batch}
  if [[ -e "${final}" ]]; then echo "batch ${batch} already final; skipping"; return 0; fi
  if [[ ! -d "${pending}" ]]; then echo "batch ${batch} pending dir absent: ${pending}" >&2; return 1; fi
  python "${PACKAGE}/22_validate_replay_batch.py" --batch "${batch}" --pending "${pending}"
  mkdir -p "$(dirname "${final}")"
  mv -- "${pending}" "${final}"
  echo "batch ${batch} promoted to ${final}"
}

revalidate 3 "${GEN}/.batch_3.pending.20918337"
revalidate 4 "${GEN}/.batch_4.pending.20918351"
echo "=== replay_execution ==="
ls "${GEN}/replay_execution"
