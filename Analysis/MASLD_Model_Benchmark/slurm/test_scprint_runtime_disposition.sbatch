#!/usr/bin/env bash
#SBATCH --job-name=model-test-305
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
cd "${root}"
export PYTHONPATH="${root}/src"
/gpfs/commons/home/jameslee/micromamba/envs/decima/bin/python -m unittest -v \
  tests.unit.test_scprint_runtime_disposition \
  tests.contract.test_registry_files.RegistryFileContractTest.test_all_registry_files_parse_and_have_expected_schema \
  tests.contract.test_registry_files.RegistryFileContractTest.test_model_roster_is_complete_unique_and_fail_closed \
  tests.contract.test_registry_files.RegistryFileContractTest.test_scprint_v1_preflight_rejects_silent_replacement
