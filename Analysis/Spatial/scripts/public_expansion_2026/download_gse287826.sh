#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
DATA_DIR="${PROJECT_ROOT}/Analysis/Spatial/data/public_expansion_2026/gse287826"
CANDIDATE_ID="program-context-v2-candidate-2026-08-07"
AUDIT_DIR="${PROJECT_ROOT}/Analysis/Spatial/candidates/${CANDIDATE_ID}/acquisition/gse287826"
MANIFEST="${AUDIT_DIR}/download_manifest.tsv"

mkdir -p "${DATA_DIR}" "${AUDIT_DIR}"

download_file() {
  local filename="$1"
  local url="$2"
  local expected_bytes="$3"
  local expected_sha="$4"
  local destination="${DATA_DIR}/${filename}"
  local partial="${destination}.part"
  local headers="${AUDIT_DIR}/${filename}.headers.txt"

  if [[ ! -f "${destination}" ]]; then
    curl --fail --location --retry 12 --retry-delay 10 --retry-all-errors \
      --continue-at - --dump-header "${headers}" --output "${partial}" "${url}"
    mv "${partial}" "${destination}"
  fi

  local observed_bytes observed_sha retrieved
  observed_bytes="$(stat -c '%s' "${destination}")"
  observed_sha="$(sha256sum "${destination}" | awk '{print $1}')"
  retrieved="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

  if [[ "${expected_bytes}" != "NA" && "${observed_bytes}" != "${expected_bytes}" ]]; then
    echo "Byte-size mismatch for ${filename}: ${observed_bytes} != ${expected_bytes}" >&2
    exit 1
  fi
  if [[ "${expected_sha}" != "NA" && "${observed_sha}" != "${expected_sha}" ]]; then
    echo "SHA256 mismatch for ${filename}: ${observed_sha} != ${expected_sha}" >&2
    exit 1
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${filename}" "${url}" "${expected_bytes}" "${observed_bytes}" \
    "${expected_sha}" "${observed_sha}" "${retrieved}" >> "${MANIFEST}"
}

printf 'file\turl\texpected_bytes\tobserved_bytes\texpected_sha256\tobserved_sha256\tretrieved_utc\n' > "${MANIFEST}"

download_file \
  "GSE287826_All_Data_WTA.xlsx" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287826/suppl/GSE287826_All_Data_WTA.xlsx" \
  "8392802" \
  "d805dbea79e6a88e7339e5fbfa44c4e7fb68b2226f69a5e570cb353b897007f1"

download_file \
  "GSE287826_family.soft.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287826/soft/GSE287826_family.soft.gz" \
  "NA" \
  "NA"

SPATIAL_PYTHON="${SPATIAL_PYTHON:-/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python}"
test -x "${SPATIAL_PYTHON}" || {
  echo "Missing spatial-environment Python: ${SPATIAL_PYTHON}" >&2
  exit 2
}
"${SPATIAL_PYTHON}" \
  "${PROJECT_ROOT}/Analysis/Spatial/scripts/public_expansion_2026/audit_gse287826_donor_gate.py" \
  --workbook "${DATA_DIR}/GSE287826_All_Data_WTA.xlsx" \
  --soft "${DATA_DIR}/GSE287826_family.soft.gz" \
  --output-dir "${AUDIT_DIR}"
