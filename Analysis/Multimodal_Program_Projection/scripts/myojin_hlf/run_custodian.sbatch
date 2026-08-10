#!/usr/bin/env bash
#SBATCH --job-name=xlsx
#SBATCH --partition=io
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/logs/xlsx_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/logs/xlsx_%j.err

set -euo pipefail
PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
CANDIDATE_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf"
SCRIPT_DIR="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/myojin_hlf"
RAW="${CANDIDATE_ROOT}/source/raw/Myojin_Supplemental_Table_S1_C218.xlsx"

mkdir -p "${CANDIDATE_ROOT}/source/raw" "${CANDIDATE_ROOT}/logs"
if [[ ! -f "${RAW}" ]]; then
  curl -L --fail --retry 5 --connect-timeout 30 --output "${RAW}.partial" https://links.lww.com/HC9/C218
  mv "${RAW}.partial" "${RAW}"
  chmod 440 "${RAW}"
fi
python3 -B "${SCRIPT_DIR}/00_mask_myojin_workbook.py" \
  --workbook "${RAW}" --outdir "${CANDIDATE_ROOT}" \
  --source-url https://links.lww.com/HC9/C218
