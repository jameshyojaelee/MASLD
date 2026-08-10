#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
CANDIDATE_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf"
SOURCE_DIR="${CANDIDATE_ROOT}/source/depmap_24q4"
SCRIPT_DIR="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/scripts/myojin_hlf"

mkdir -p "${SOURCE_DIR}"

download_verified() {
  local file_id="$1"
  local filename="$2"
  local expected_md5="$3"
  local destination="${SOURCE_DIR}/${filename}"
  local partial="${destination}.partial"

  if [[ -f "${destination}" ]] && [[ "$(md5sum "${destination}" | awk '{print $1}')" == "${expected_md5}" ]]; then
    return 0
  fi
  curl -L --fail --retry 5 --retry-delay 10 --connect-timeout 30 \
    --output "${partial}" "https://ndownloader.figshare.com/files/${file_id}"
  local observed_md5
  observed_md5="$(md5sum "${partial}" | awk '{print $1}')"
  if [[ "${observed_md5}" != "${expected_md5}" ]]; then
    echo "MD5 mismatch for ${filename}: ${observed_md5} != ${expected_md5}" >&2
    exit 1
  fi
  mv "${partial}" "${destination}"
  chmod 440 "${destination}"
}

download_verified 51065297 Model.csv 675210d17675f3517b0ce39a3c274f16
download_verified 51064667 CRISPRGeneEffect.csv 6edf7ade09b9b34199210b559d4745d3
download_verified 51065489 OmicsExpressionProteinCodingGenesTPMLogp1.csv 71794802b750ce77c422dad0720a40af

python3 -B "${SCRIPT_DIR}/01_extract_depmap_hlf.py" \
  --source-dir "${SOURCE_DIR}" \
  --outdir "${CANDIDATE_ROOT}" \
  --depmap-release "DepMap Public 24Q4"
