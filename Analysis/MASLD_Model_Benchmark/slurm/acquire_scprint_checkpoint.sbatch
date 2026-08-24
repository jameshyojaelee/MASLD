#!/usr/bin/env bash
#SBATCH --job-name=masld-scprint-acquire
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

package_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
target="${package_root}/executions/scprint-acquisition-${SLURM_JOB_ID}"
stage="${target}.staging"
checkpoint="${stage}/v2-medium.ckpt"
partial="${checkpoint}.part"
expected_bytes="221775592"
expected_sha256="a4cf0753270d4ff451a5dbddadd4c59a53e10e4c3e96a8af0db407ad893c36c5"
source_url="https://huggingface.co/jkobject/scPRINT/resolve/f60c87152e83fedb932a84f8ba850fbf489a5a86/v2-medium.ckpt"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite scPRINT acquisition: %s\n' "${target}" >&2
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
  "${source_url}"

observed_bytes="$(stat --format='%s' "${partial}")"
if [[ "${observed_bytes}" != "${expected_bytes}" ]]; then
  printf 'Byte-size mismatch: expected %s, observed %s\n' \
    "${expected_bytes}" "${observed_bytes}" >&2
  exit 1
fi
printf '%s  %s\n' "${expected_sha256}" "${partial}" | sha256sum -c -
mv "${partial}" "${checkpoint}"
(
  cd "${stage}"
  sha256sum v2-medium.ckpt > SHA256SUMS
)
chmod 0440 "${checkpoint}" "${stage}/SHA256SUMS"
mv "${stage}" "${target}"
printf 'scPRINT acquisition candidate: %s\n' "${target}"
printf 'No checkpoint content was deserialized or executed. Restricted inspection remains mandatory.\n'
