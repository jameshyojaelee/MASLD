#!/usr/bin/env bash
#SBATCH --job-name=corgi-plus-disposition
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "${ROOT}"

for pair in \
  "43210eabff99797a7c026e52f5ec19ead336424ba0ff13a237e71203232a173f config/artifacts/models/corgi_plus/release_disposition_20260824.json" \
  "a2fe2db2c4b7a45c67d65674db8e69a001b08cd582e47cc98c0680acddc038ac scripts/audit_corgi_plus_disposition.py" \
  "3c023b2cbf137626fc1366f228bef9e5a290920b6d0ca68baaafccf2b7c9b7b5 tests/unit/test_audit_corgi_plus_disposition.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest \
  tests.unit.test_audit_corgi_plus_disposition
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_corgi_plus_disposition.py
