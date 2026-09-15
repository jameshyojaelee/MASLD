#!/usr/bin/env bash
#SBATCH --job-name=remaining-multimodal-audit
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
  "de07b243f55b38a616fee8c333540a50baa2ec284f2d08944b6f5a259171e591 config/artifacts/models/remaining_multimodal_families/disposition_20260824.json" \
  "2bb53407ad32b4247f9e7a39d9a05028f282a903b1e8d57960b0478aec00bcaf scripts/audit_remaining_multimodal_families_disposition.py" \
  "a51a62e4c69157e542a08d3fec813aed772299cfed17c365a725a056942d5678 tests/unit/test_audit_remaining_multimodal_families_disposition.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest \
  tests.unit.test_audit_remaining_multimodal_families_disposition
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_remaining_multimodal_families_disposition.py
