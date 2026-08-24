#!/usr/bin/env bash
#SBATCH --job-name=masld-scgpt-acquire
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

package_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
target="${package_root}/executions/scgpt-acquisition-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite scGPT acquisition: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 \
  "${stage}/scgpt_whole_human" \
  "${stage}/scgpt_continual"

download_drive_file() {
  local file_id="$1"
  local destination="$2"
  local expected_bytes="$3"
  local expected_sha256="$4"
  local partial="${destination}.part"
  if [[ -e "${destination}" || -e "${partial}" ]]; then
    printf 'Refusing to overwrite downloaded artifact: %s\n' "${destination}" >&2
    exit 1
  fi
  curl \
    --fail \
    --location \
    --retry 5 \
    --retry-all-errors \
    --show-error \
    --silent \
    --output "${partial}" \
    "https://drive.usercontent.google.com/download?id=${file_id}&export=download&confirm=t"
  local observed_bytes
  observed_bytes="$(stat --format='%s' "${partial}")"
  if [[ "${observed_bytes}" != "${expected_bytes}" ]]; then
    printf 'Byte-size mismatch for %s: expected %s, observed %s\n' \
      "${destination}" "${expected_bytes}" "${observed_bytes}" >&2
    exit 1
  fi
  if [[ "${expected_sha256}" != "UNRESOLVED" ]]; then
    printf '%s  %s\n' "${expected_sha256}" "${partial}" | sha256sum -c -
  fi
  mv "${partial}" "${destination}"
}

download_drive_file \
  "14AebJfGOUF047Eg40hk57HCtrb0fyDTm" \
  "${stage}/scgpt_whole_human/best_model.pt" \
  "205385258" \
  "UNRESOLVED"
download_drive_file \
  "1hh2zGKyWAx3DyovD30GStZ3QlzmSqdk1" \
  "${stage}/scgpt_whole_human/args.json" \
  "1300" \
  "c18e075e018140cb8b2d9029387b9de26607a5ce6a8ccabd6ead70cd76b95d60"
download_drive_file \
  "1H3E_MJ-Dl36AQV6jLbna2EdvgPaqvqcC" \
  "${stage}/scgpt_whole_human/vocab.json" \
  "1317639" \
  "acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f"
download_drive_file \
  "1x1SfmFdI-zcocmqWAd7ZTC9CTEAVfKZq" \
  "${stage}/scgpt_continual/best_model.pt" \
  "207861754" \
  "UNRESOLVED"
download_drive_file \
  "15TEZmd2cZCrHwgfE424fgQkGUZCXiYrR" \
  "${stage}/scgpt_continual/args.json" \
  "1503" \
  "77fec83c32306225f37cc023321b716a49ce9e778039f2c1c4b80e1e4e7008bd"
download_drive_file \
  "1jfT_T5n8WNbO9QZcLWObLdRG8lYFKH-Q" \
  "${stage}/scgpt_continual/vocab.json" \
  "1317639" \
  "ee2b2c90158eedb97c2318e49abaaed0a02c6fdf7e3f7ca6a821906413c4d2a4"

(
  cd "${stage}"
  find scgpt_whole_human scgpt_continual -type f -print0 \
    | sort -z \
    | xargs -0 sha256sum > SHA256SUMS
)
chmod 0440 "${stage}/SHA256SUMS"
find "${stage}/scgpt_whole_human" "${stage}/scgpt_continual" \
  -type f -exec chmod 0440 {} +
mv "${stage}" "${target}"
printf 'scGPT acquisition candidate: %s\n' "${target}"
printf 'No checkpoint content was deserialized or executed.\n'
