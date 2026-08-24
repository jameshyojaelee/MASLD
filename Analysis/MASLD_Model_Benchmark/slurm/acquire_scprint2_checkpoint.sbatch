#!/usr/bin/env bash
#SBATCH --job-name=masld-scprint2-acquire
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

package_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
target="${package_root}/executions/scprint2-acquisition-${SLURM_JOB_ID}"
stage="${target}.staging"
checkpoint="${stage}/18hebyht-final-small.ckpt"
partial="${checkpoint}.part"
expected_bytes="889706878"
expected_sha256="2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822"
source_url="https://huggingface.co/jkobject/scPRINT/resolve/d02d4aa684a8422d8bfdc1dbb521a85c773a2772/18hebyht-final-small.ckpt"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite scPRINT-2 acquisition: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"
curl --fail --location --retry 5 --retry-all-errors --show-error --silent \
  --output "${partial}" "${source_url}"
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
  sha256sum 18hebyht-final-small.ckpt > SHA256SUMS
)
chmod 0440 "${checkpoint}" "${stage}/SHA256SUMS"
mv "${stage}" "${target}"
printf 'scPRINT-2 acquisition candidate: %s\n' "${target}"
printf 'No checkpoint content was deserialized or executed. Restricted inspection remains mandatory.\n'
