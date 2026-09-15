#!/usr/bin/env bash
#SBATCH --job-name=model-test-305
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
python="/gpfs/commons/home/jameslee/micromamba/envs/gca-test/bin/python"
export PYTHONPATH="${root}/src"
export PYTHONDONTWRITEBYTECODE=1
cd "${root}"

"${python}" -m unittest -v \
  tests.contract.test_registry_files.RegistryFileContractTest.test_model_roster_is_complete_unique_and_fail_closed \
  tests.contract.test_registry_files.RegistryFileContractTest.test_scprint2_preflight_fails_closed_on_exposure_and_neighbors \
  tests.unit.test_scprint2_embedding_contract
