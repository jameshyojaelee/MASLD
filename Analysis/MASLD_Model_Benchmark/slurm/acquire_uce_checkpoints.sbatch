#!/usr/bin/env bash
#SBATCH --job-name=masld-uce-acquire
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

package_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
target="${package_root}/executions/uce-acquisition-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite UCE acquisition: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

download_figshare_file() {
  local file_id="$1"
  local destination="$2"
  local expected_bytes="$3"
  local expected_md5="$4"
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
    "https://ndownloader.figshare.com/files/${file_id}"
  local observed_bytes
  observed_bytes="$(stat --format='%s' "${partial}")"
  if [[ "${observed_bytes}" != "${expected_bytes}" ]]; then
    printf 'Byte-size mismatch for %s: expected %s, observed %s\n' \
      "${destination}" "${expected_bytes}" "${observed_bytes}" >&2
    exit 1
  fi
  printf '%s  %s\n' "${expected_md5}" "${partial}" | md5sum -c -
  mv "${partial}" "${destination}"
}

download_figshare_file 42706576 "${stage}/4layer_model.torch" 3403514339 3e2f59d6da6eaa5396aa297edfcec08f
download_figshare_file 43423236 "${stage}/33l_8ep_1024t_1280.torch" 5686283829 eed62d4fe51873c87bcf9660faa07707
download_figshare_file 42706555 "${stage}/species_offsets.pkl" 139 ddf9143508857e62607f323b338ff383
download_figshare_file 42706558 "${stage}/species_chrom.csv" 4097781 42a9b5871b1ed955fd05414b3ab89081
download_figshare_file 42706585 "${stage}/all_tokens.torch" 2979205876 3523ea431c403f95701b6eeb0c30b997
download_figshare_file 42715213 "${stage}/protein_embeddings.tar.gz" 2735410523 6756d7fa8348e0f70b5e5b11e324f023

(
  cd "${stage}"
  sha256sum \
    4layer_model.torch \
    33l_8ep_1024t_1280.torch \
    species_offsets.pkl \
    species_chrom.csv \
    all_tokens.torch \
    protein_embeddings.tar.gz \
    > SHA256SUMS
)
chmod 0440 "${stage}/SHA256SUMS"
find "${stage}" -type f -exec chmod 0440 {} +
mv "${stage}" "${target}"
printf 'UCE acquisition candidate: %s\n' "${target}"
printf 'No checkpoint content was deserialized, executed, or extracted.\n'
