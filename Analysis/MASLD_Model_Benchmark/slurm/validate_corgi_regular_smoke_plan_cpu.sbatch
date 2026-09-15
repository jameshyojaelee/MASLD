#!/usr/bin/env bash
#SBATCH --job-name=corgi-regular-smoke-audit
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
  "92543e0ee28e0cdbdd9e2bb13a32b0b4f2960a9320dc417aa881d4d1f8f63635 config/artifacts/models/corgi/bounded_prediction_smoke_disposition_20260824.json" \
  "ff2c07ca2bd5c1d38da7c36682826c3fa7b7fa1e6fc76d6bcfc648c16376b668 scripts/audit_corgi_regular_smoke_plan.py" \
  "1d15bf6d1a06a68921aa9bbecd4be744276ca4b8feaff12073995eee1ffcb64b tests/unit/test_audit_corgi_regular_smoke_plan.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest \
  tests.unit.test_audit_corgi_regular_smoke_plan
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_corgi_regular_smoke_plan.py
