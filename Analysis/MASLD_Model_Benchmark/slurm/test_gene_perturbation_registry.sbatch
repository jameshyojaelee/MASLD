#!/usr/bin/env bash
#SBATCH --job-name=masld-gene-perturb-registry
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/gene-perturbation-registry-${SLURM_JOB_ID}"
if [[ -e "${TARGET}" ]]; then
  printf 'Refusing to overwrite registry test: %s\n' "${TARGET}" >&2
  exit 1
fi
mkdir "${TARGET}"
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python -m unittest -v \
  tests.contract.test_registry_files.RegistryFileContractTest.test_all_registry_files_parse_and_have_expected_schema \
  tests.contract.test_registry_files.RegistryFileContractTest.test_model_roster_is_complete_unique_and_fail_closed \
  tests.contract.test_registry_files.RegistryFileContractTest.test_every_model_audit_bundle_is_complete_and_hash_bound \
  tests.contract.test_registry_files.RegistryFileContractTest.test_perturbation_models_reject_incompatible_held_targets \
  tests.contract.test_registry_files.RegistryFileContractTest.test_remaining_perturbation_models_reject_unseen_noncoding_targets \
  tests.contract.test_registry_files.RegistryFileContractTest.test_task_dataset_and_split_references_resolve \
  tests.contract.test_registry_files.RegistryFileContractTest.test_split_contract_uses_biological_units \
  2>&1 | tee "${TARGET}/unittest.log"
PYTHONPATH=src python - <<'PY' | tee "${TARGET}/strict_registry.log"
from pathlib import Path
from masld_bench.registry import Registry

registry = Registry.load(Path("config"), validate_references=False)
print(len(registry.datasets), len(registry.models), len(registry.tasks), len(registry.splits))
print(registry.tasks["gene_perturbation_response"].to_json())
assert registry.tasks["gene_perturbation_response"].primary_evaluator_id == "not_applicable_nonpromotable_v1"
assert registry.tasks["gene_perturbation_response"].datasets_sealed == ()
PY
sha256sum \
  config/datasets/gse238219.toml \
  config/datasets/gse264667.toml \
  config/datasets/gse242934.toml \
  config/datasets/cra009621_oscar.toml \
  config/splits/gene_perturbation_outer.toml \
  config/tasks/gene_perturbation_response.toml \
  config/models/perturbation.toml \
  config/models/mandatory_baselines.toml \
  src/masld_bench/registry.py \
  tests/contract/test_registry_files.py \
  > "${TARGET}/source.sha256"
PYTHONPATH=src python - "${TARGET}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

print(freeze_tree(Path(sys.argv[1]), {
    "artifact_class": "gene_perturbation_registry_validation",
    "gse281160_remains_noncoding_task_only": True,
    "gse313774_accessed": False,
    "biological_replicate_policy": "cells_guides_lanes_files_not_replicates",
}))
PY
sha256sum "${TARGET}/ARTIFACTS.json"
