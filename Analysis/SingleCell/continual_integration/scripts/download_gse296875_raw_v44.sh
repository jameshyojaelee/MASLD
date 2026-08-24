#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: download_gse296875_raw_v44.sh SOURCE_TSV OUTPUT_DIRECTORY" >&2
  exit 2
fi

source_tsv=$(realpath "$1")
output_directory=$(realpath -m "$2")
mkdir -p -- "$output_directory"

download_one() {
  local well=$1
  local gsm=$2
  local url=$3
  local final="${GSE296875_OUTPUT_DIRECTORY}/${gsm}_${well}_raw_feature_bc_matrix.h5"
  local partial="${final}.partial"
  if [[ -e "$final" ]]; then
    echo "refusing to overwrite ${final}" >&2
    return 1
  fi
  curl --location --fail --silent --show-error --retry 4 --retry-delay 5 --continue-at - \
    --output "$partial" "$url"
  mv -- "$partial" "$final"
  sha256sum "$final" > "${final}.sha256"
}
export -f download_one
export GSE296875_OUTPUT_DIRECTORY="$output_directory"

tail -n +2 "$source_tsv" | xargs -P 4 -n 3 bash -c \
  'download_one "$1" "$2" "$3"' _
