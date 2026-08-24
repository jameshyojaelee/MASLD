#!/usr/bin/env bash
#SBATCH --job-name=masld-scim-acquire
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

package_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
target="${package_root}/executions/scimilarity-acquisition-${SLURM_JOB_ID}"
stage="${target}.staging"
archive="${stage}/model_v1.1.tar.gz"
partial="${archive}.part"
expected_bytes="30310810843"
expected_md5="546251b7c435f3b1dbe38e2e420ad57f"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite SCimilarity acquisition: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

curl \
  --fail \
  --location \
  --retry 5 \
  --retry-all-errors \
  --show-error \
  --silent \
  --output "${partial}" \
  "https://zenodo.org/api/records/10685499/files/model_v1.1.tar.gz/content"

observed_bytes="$(stat --format='%s' "${partial}")"
if [[ "${observed_bytes}" != "${expected_bytes}" ]]; then
  printf 'Byte-size mismatch: expected %s, observed %s\n' \
    "${expected_bytes}" "${observed_bytes}" >&2
  exit 1
fi
printf '%s  %s\n' "${expected_md5}" "${partial}" | md5sum -c -
mv "${partial}" "${archive}"
(
  cd "${stage}"
  sha256sum model_v1.1.tar.gz > SHA256SUMS
)
chmod 0440 "${archive}" "${stage}/SHA256SUMS"
mv "${stage}" "${target}"
printf 'SCimilarity acquisition candidate: %s\n' "${target}"
printf 'The archive was not extracted and checkpoint content was not deserialized or executed.\n'
