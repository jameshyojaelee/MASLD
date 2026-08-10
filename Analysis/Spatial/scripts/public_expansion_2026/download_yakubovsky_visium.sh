#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
DATA_DIR="${PROJECT_ROOT}/Analysis/Spatial/data/public_expansion_2026/yakubovsky2026"
CANDIDATE_ID="program-context-v2-candidate-2026-08-07"
AUDIT_DIR="${PROJECT_ROOT}/Analysis/Spatial/candidates/${CANDIDATE_ID}/acquisition/yakubovsky2026"
DESTINATION="${DATA_DIR}/Visium.zip"
PARTIAL="${DESTINATION}.part"
HEADERS="${AUDIT_DIR}/Visium.zip.headers.txt"
MANIFEST="${AUDIT_DIR}/visium_download_manifest.tsv"
URL="https://zenodo.org/api/records/17735506/files/Visium.zip/content"
EXPECTED_BYTES="30479914954"
EXPECTED_MD5="8e31a754050ba82ace66246dfb661201"

mkdir -p "${DATA_DIR}" "${AUDIT_DIR}"

if [[ ! -f "${DESTINATION}" ]]; then
  curl --fail --location --retry 12 --retry-delay 10 --retry-all-errors \
    --continue-at - --dump-header "${HEADERS}" --output "${PARTIAL}" "${URL}"
  mv "${PARTIAL}" "${DESTINATION}"
fi

OBSERVED_BYTES="$(stat -c '%s' "${DESTINATION}")"
OBSERVED_MD5="$(md5sum "${DESTINATION}" | awk '{print $1}')"
OBSERVED_SHA="$(sha256sum "${DESTINATION}" | awk '{print $1}')"
RETRIEVED="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

[[ "${OBSERVED_BYTES}" == "${EXPECTED_BYTES}" ]] || {
  echo "Byte-size mismatch for Visium.zip: ${OBSERVED_BYTES} != ${EXPECTED_BYTES}" >&2
  exit 1
}
[[ "${OBSERVED_MD5}" == "${EXPECTED_MD5}" ]] || {
  echo "MD5 mismatch for Visium.zip: ${OBSERVED_MD5} != ${EXPECTED_MD5}" >&2
  exit 1
}

printf 'file\turl\texpected_bytes\tobserved_bytes\texpected_md5\tobserved_md5\tsha256\tretrieved_utc\textraction_status\n' > "${MANIFEST}"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "Visium.zip" "${URL}" "${EXPECTED_BYTES}" "${OBSERVED_BYTES}" \
  "${EXPECTED_MD5}" "${OBSERVED_MD5}" "${OBSERVED_SHA}" "${RETRIEVED}" \
  "not_extracted_pending_lean_schema_gate" >> "${MANIFEST}"
