#!/usr/bin/env bash
#SBATCH --job-name=proteomics
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Proteomics/logs/integration_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Proteomics/logs/integration_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT="${PROJECT_ROOT}/Analysis/Proteomics/scripts/proteo_transcriptomic_integration.py"

mkdir -p "${PROJECT_ROOT}/Analysis/Proteomics/logs"

eval "$(micromamba shell hook -s bash)"
set +u  # conda activation references unbound vars (ADDR2LINE) — aborts under set -u
micromamba activate rnaseq
set -u

export PYTHONNOUSERSITE=1

python3 "${SCRIPT}" --project-root "${PROJECT_ROOT}"
