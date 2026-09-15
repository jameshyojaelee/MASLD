#!/usr/bin/env bash
#SBATCH --job-name=observed-multiome-audit
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
  "a474db7f50fe481206b017c5d77b1434d39f48c55683b3e2a4c9a544229ec0b3 config/artifacts/models/observed_multiome_integration/disposition_20260824.json" \
  "8bb1027133b4216481c5b61fa29b266a394ef13b1dea9d74b1c6667ce66e78ec scripts/audit_observed_multiome_integration_disposition.py" \
  "d963e4d60bc71cb2cf4979349f38addd097e2953b19a837f566f1f72691256a2 tests/unit/test_audit_observed_multiome_integration_disposition.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest \
  tests.unit.test_audit_observed_multiome_integration_disposition
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_observed_multiome_integration_disposition.py
