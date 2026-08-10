#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
DATA_DIR="${PROJECT_ROOT}/Analysis/Spatial/data/public_expansion_2026/yakubovsky2026"
CANDIDATE_ID="program-context-v2-candidate-2026-08-07"
AUDIT_DIR="${PROJECT_ROOT}/Analysis/Spatial/candidates/${CANDIDATE_ID}/acquisition/yakubovsky2026"
MANIFEST="${AUDIT_DIR}/lean_download_manifest.tsv"

mkdir -p "${DATA_DIR}" "${AUDIT_DIR}"
printf 'file\turl\texpected_bytes\tobserved_bytes\texpected_md5\tobserved_md5\tsha256\tretrieved_utc\n' > "${MANIFEST}"

download_md5() {
  local filename="$1"
  local url="$2"
  local expected_bytes="$3"
  local expected_md5="$4"
  local destination="${DATA_DIR}/${filename}"
  local partial="${destination}.part"
  local headers="${AUDIT_DIR}/${filename}.headers.txt"

  if [[ ! -f "${destination}" ]]; then
    curl --fail --location --retry 12 --retry-delay 10 --retry-all-errors \
      --continue-at - --dump-header "${headers}" --output "${partial}" "${url}"
    mv "${partial}" "${destination}"
  fi

  local observed_bytes observed_md5 observed_sha retrieved
  observed_bytes="$(stat -c '%s' "${destination}")"
  observed_md5="$(md5sum "${destination}" | awk '{print $1}')"
  observed_sha="$(sha256sum "${destination}" | awk '{print $1}')"
  retrieved="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

  [[ "${observed_bytes}" == "${expected_bytes}" ]] || {
    echo "Byte-size mismatch for ${filename}: ${observed_bytes} != ${expected_bytes}" >&2
    exit 1
  }
  [[ "${observed_md5}" == "${expected_md5}" ]] || {
    echo "MD5 mismatch for ${filename}: ${observed_md5} != ${expected_md5}" >&2
    exit 1
  }

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${filename}" "${url}" "${expected_bytes}" "${observed_bytes}" \
    "${expected_md5}" "${observed_md5}" "${observed_sha}" "${retrieved}" >> "${MANIFEST}"
}

download_md5 \
  "human_samples_metadata.xlsx" \
  "https://zenodo.org/api/records/17735506/files/human_samples_metadata.xlsx/content" \
  "21271" \
  "ff13c3f429d1a3f84a6035a03e72639f"

download_md5 \
  "v.mat" \
  "https://zenodo.org/api/records/17735587/files/v.mat/content" \
  "3996202856" \
  "73f2ae74d2984363511d063af2873b0a"

download_md5 \
  "zon_struct_all_full.mat" \
  "https://zenodo.org/api/records/17735587/files/zon_struct_all_full.mat/content" \
  "203729828" \
  "8dc2e38d58c84e16a5143c9a60d02146"

CODE_DIR="${DATA_DIR}/Human-liver"
if [[ ! -d "${CODE_DIR}/.git" ]]; then
  git clone --depth 1 https://github.com/OranYak/Human-liver.git "${CODE_DIR}"
fi
git -C "${CODE_DIR}" rev-parse HEAD > "${AUDIT_DIR}/github_commit.txt"
git -C "${CODE_DIR}" remote get-url origin > "${AUDIT_DIR}/github_origin.txt"

SPATIAL_PYTHON="${SPATIAL_PYTHON:-/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python}"
test -x "${SPATIAL_PYTHON}" || {
  echo "Missing spatial-environment Python: ${SPATIAL_PYTHON}" >&2
  exit 2
}
"${SPATIAL_PYTHON}" \
  "${PROJECT_ROOT}/Analysis/Spatial/scripts/public_expansion_2026/audit_yakubovsky_schema.py" \
  --data-dir "${DATA_DIR}" \
  --output-dir "${AUDIT_DIR}"
