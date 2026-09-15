#!/usr/bin/env bash
#SBATCH --job-name=masld-tf-acquire
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

package_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
target="${package_root}/executions/transcriptformer-acquisition-${SLURM_JOB_ID}"
stage="${target}.staging"
archive="${stage}/tf_sapiens.tar.gz"
partial="${archive}.part"
response_headers="${stage}/response.headers"
version_id="P5N5hx41EeWC7I0l8oW7yDyiJHS5InPU"
expected_etag='"a841d3efcfc21eb14c41c6151bc39b73-217"'
expected_bytes="1819975447"
source_url="https://czi-transcriptformer.s3.amazonaws.com/weights/tf_sapiens.tar.gz?versionId=${version_id}"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer acquisition: %s\n' "${target}" >&2
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
  --dump-header "${response_headers}" \
  --output "${partial}" \
  "${source_url}"

observed_version="$(tr -d '\r' < "${response_headers}" | \
  sed -n 's/^[Xx]-[Aa][Mm][Zz]-[Vv][Ee][Rr][Ss][Ii][Oo][Nn]-[Ii][Dd]:[[:space:]]*//p' | \
  tail -n 1)"
observed_etag="$(tr -d '\r' < "${response_headers}" | \
  sed -n 's/^[Ee][Tt][Aa][Gg]:[[:space:]]*//p' | tail -n 1)"
if [[ "${observed_version}" != "${version_id}" ]]; then
  printf 'S3 version mismatch: expected %s, observed %s\n' \
    "${version_id}" "${observed_version}" >&2
  exit 1
fi
if [[ "${observed_etag}" != "${expected_etag}" ]]; then
  printf 'S3 ETag mismatch: expected %s, observed %s\n' \
    "${expected_etag}" "${observed_etag}" >&2
  exit 1
fi

observed_bytes="$(stat --format='%s' "${partial}")"
if [[ "${observed_bytes}" != "${expected_bytes}" ]]; then
  printf 'Byte-size mismatch: expected %s, observed %s\n' \
    "${expected_bytes}" "${observed_bytes}" >&2
  exit 1
fi
mv "${partial}" "${archive}"
(
  cd "${stage}"
  sha256sum tf_sapiens.tar.gz response.headers > SHA256SUMS
)
chmod 0440 "${archive}" "${response_headers}" "${stage}/SHA256SUMS"
mv "${stage}" "${target}"
printf 'TranscriptFormer acquisition candidate: %s\n' "${target}"
printf 'The archive was not extracted and checkpoint content was not deserialized or executed.\n'
