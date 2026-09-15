#!/usr/bin/env bash
#SBATCH --job-name=source-contract-256398
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "${ROOT}"
module purge
module load h5py/3.9.0-foss-2023a
PYTHONPATH=${ROOT}:${ROOT}/src:${PYTHONPATH:-} python -m unittest -v \
  tests.unit.test_audit_gse256398_h5 \
  tests.unit.test_build_gse256398_gencode49_crosswalk \
  tests.contract.test_cross_cohort_plan \
  tests.contract.test_registry_files.RegistryFileContractTest.test_all_registry_files_parse_and_have_expected_schema \
  tests.contract.test_registry_files.RegistryFileContractTest.test_task_dataset_and_split_references_resolve \
  tests.contract.test_registry_files.RegistryFileContractTest.test_pairing_topology_prevents_false_pairs \
  tests.unit.test_firewall.DatasetAccessTests.test_access_is_derived_from_a_strict_manifest
PYTHONPATH=${ROOT}/src python - <<'PY'
from pathlib import Path
from masld_bench.artifacts import verify_frozen_tree

root = Path.cwd()
for relative in (
    "executions/gse256398-source-h5-21069407",
    "executions/gse256398-gencode49-crosswalk-21069439",
):
    verify_frozen_tree(root / relative)
print("GSE256398 frozen source and crosswalk verification passed")
PY
printf '%s  %s\n' \
  2aec756d6ed95e31a16de9191f3c5ba7bfb44e3321b70661a8f7eb4bc9251fdb \
  config/tasks/cell_state_mapping.toml \
  | sha256sum --check --strict
sha256sum \
  config/datasets/gse256398.toml \
  config/evaluation/cross_cohort_expansion.toml \
  config/tasks/cell_state_mapping.toml \
  scripts/audit_gse256398_h5.py \
  scripts/build_gse256398_gencode49_crosswalk.py \
  tests/unit/test_audit_gse256398_h5.py \
  tests/unit/test_build_gse256398_gencode49_crosswalk.py \
  tests/contract/test_cross_cohort_plan.py \
  tests/contract/test_registry_files.py \
  slurm/validate_gse256398_activation.sbatch
