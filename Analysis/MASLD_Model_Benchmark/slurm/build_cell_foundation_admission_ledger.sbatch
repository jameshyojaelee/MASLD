#!/usr/bin/env bash
#SBATCH --job-name=masld-cellfm-ledger
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONNOUSERSITE=1

root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
target="${root}/executions/cell-foundation-admission-ledger-${SLURM_JOB_ID}"
stage="${target}.staging"
if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite admission ledger: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

python3 "${root}/scripts/build_cell_foundation_admission_ledger.py" \
  --root "${root}" --output "${stage}/admission_ledger.json"
python3 --version > "${stage}/python_version.txt" 2>&1
python3 -m pip freeze > "${stage}/pip_freeze.txt"
sha256sum \
  "${root}/scripts/build_cell_foundation_admission_ledger.py" \
  "${root}/slurm/build_cell_foundation_admission_ledger.sbatch" \
  > "${stage}/source_manifest.sha256"
(
  cd "${stage}"
  sha256sum admission_ledger.json python_version.txt pip_freeze.txt > SHA256SUMS
)
chmod 0440 "${stage}"/*
mv "${stage}" "${target}"
printf 'Cell-foundation admission ledger: %s\n' "${target}"
