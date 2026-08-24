#!/usr/bin/env bash
#SBATCH --job-name=masld-tf-inventory
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
archive="${root}/executions/transcriptformer-acquisition-21066078/tf_sapiens.tar.gz"
expected_sha256="d108a7a9ee7fdf38e0dcd4e58593e94304cba927df6ec63e82ea74e2b177e529"
expected_size="1819975447"
target="${root}/executions/transcriptformer-safe-inventory-${SLURM_JOB_ID}"
stage="${target}.staging"

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite TranscriptFormer inventory: %s\n' "${target}" >&2
  exit 1
fi
observed_sha256="$(sha256sum "${archive}" | cut -d ' ' -f 1)"
observed_size="$(stat --format='%s' "${archive}")"
if [[ "${observed_sha256}" != "${expected_sha256}" || \
      "${observed_size}" != "${expected_size}" ]]; then
  printf 'TranscriptFormer archive bytes differ\n' >&2
  exit 1
fi
install -d -m 0750 "${stage}"

python3 "${root}/scripts/epibert_safe_archive_inventory.py" \
  --archive "${archive}" \
  --output "${stage}/safe_archive_inventory.json" \
  --expected-size "${expected_size}"

python3 --version > "${stage}/python_version.txt" 2>&1
python3 -m pip freeze > "${stage}/pip_freeze.txt"
(
  cd "${stage}"
  sha256sum safe_archive_inventory.json python_version.txt pip_freeze.txt > SHA256SUMS
)
sha256sum \
  "${root}/scripts/epibert_safe_archive_inventory.py" \
  "${root}/slurm/transcriptformer_safe_archive_inventory.sbatch" \
  > "${stage}/source_manifest.sha256"
chmod 0440 "${stage}"/*
mv "${stage}" "${target}"
printf 'TranscriptFormer safe archive inventory: %s\n' "${target}"
