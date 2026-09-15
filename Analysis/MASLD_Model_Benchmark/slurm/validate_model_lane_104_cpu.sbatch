#!/usr/bin/env bash
#SBATCH --job-name=model-validate-104
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
module load Python/3.11.3-GCCcore-12.3.0

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
export PYTHONPATH="${ROOT}/src"
cd "${ROOT}"
python3 -m unittest \
  tests.unit.test_audit_epibert_epcotv2_lane \
  tests.contract.test_registry_files.RegistryFileContractTest.test_all_registry_files_parse_and_have_expected_schema \
  tests.contract.test_registry_files.RegistryFileContractTest.test_model_roster_is_complete_unique_and_fail_closed \
  tests.contract.test_registry_files.RegistryFileContractTest.test_epibert_and_epcotv2_require_observed_atac_context
