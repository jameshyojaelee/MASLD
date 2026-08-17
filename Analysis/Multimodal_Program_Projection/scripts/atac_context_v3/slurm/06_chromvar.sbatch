#!/bin/bash
#SBATCH --job-name=atacmotif
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_chromvar_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_chromvar_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
CANDIDATE=${ROOT}/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1
PREFIX=${ROOT}/.mamba/atac-context-v3-chromvar-r46
if [[ ! -x "${PREFIX}/bin/Rscript" || ! -f "${PREFIX}/chromvar.explicit.txt" ]]; then
  echo "Pinned chromVAR environment is incomplete" >&2
  exit 1
fi
if [[ -e "${CANDIDATE}/environments/chromvar.explicit.txt" ]]; then
  if ! cmp -s "${PREFIX}/chromvar.explicit.txt" "${CANDIDATE}/environments/chromvar.explicit.txt"; then
    echo "Candidate chromVAR lock differs from the pinned environment" >&2
    exit 1
  fi
else
  cp "${PREFIX}/chromvar.explicit.txt" "${CANDIDATE}/environments/chromvar.explicit.txt"
fi
"${PREFIX}/bin/Rscript" "${PACKAGE}/06_run_chromvar.R"
