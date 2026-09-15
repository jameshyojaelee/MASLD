#!/usr/bin/env bash
#SBATCH --job-name=model-validate-201
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
  tests.unit.test_audit_scooby_released_lane \
  tests.unit.scooby_test_cross_review \
  tests.unit.scooby_test_remote_safetensors_header \
  tests.unit.scooby_test_remote_zip_inventory \
  tests.contract.test_registry_files.RegistryFileContractTest.test_scooby_context_and_output_boundaries_are_frozen
