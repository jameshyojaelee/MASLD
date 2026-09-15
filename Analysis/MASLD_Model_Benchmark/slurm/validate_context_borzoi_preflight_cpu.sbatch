#!/usr/bin/env bash
#SBATCH --job-name=context-borzoi-preflight
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
  "5ed242cb79f98496a0cbffc5d5be9d34a02ba43eb3a383cd45848443afd2791d config/artifacts/models/context_borzoi/outcome_blind_architecture_preflight_20260824.json" \
  "85c78f7b10b4e23980a5b4a0a29ec54e2fab85e76445d50a9e302f2c7633114d scripts/audit_context_borzoi_preflight.py" \
  "8de6ee272b7da08a008499fb3c23e18a5692d170b410e55fe5cf34266a33ff01 tests/unit/test_audit_context_borzoi_preflight.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src:${ROOT}" python -m unittest \
  tests.unit.test_audit_context_borzoi_preflight
PYTHONPATH="${ROOT}/src:${ROOT}" python \
  scripts/audit_context_borzoi_preflight.py
