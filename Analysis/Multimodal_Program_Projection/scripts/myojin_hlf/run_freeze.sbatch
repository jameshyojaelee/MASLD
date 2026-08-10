#!/usr/bin/env bash
#SBATCH --job-name=python
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/logs/python_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/logs/python_%j.err

set -euo pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
CANDIDATE_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf"
SCRIPT_DIR="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/myojin_hlf"

python3 -B "${SCRIPT_DIR}/02_freeze_blind_spec.py" --outdir "${CANDIDATE_ROOT}"
python3 -B "${SCRIPT_DIR}/03_validate_blind_spec.py" --outdir "${CANDIDATE_ROOT}"
