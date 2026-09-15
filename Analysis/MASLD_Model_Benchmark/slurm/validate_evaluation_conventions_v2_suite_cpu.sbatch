#!/usr/bin/env bash
#SBATCH --job-name=model-check-903
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
# gpfs, not /scratch: /scratch is node-local and unreadable from the login node.
OUT=${ROOT}/executions/model-check-903-${SLURM_JOB_ID}-suite
install -d -m 0750 "${OUT}"

cd "${ROOT}"
hostname; "${PYTHON}" -V

# EXACTLY the scope PREEXISTING_FAILURES was calibrated against, plus the two
# modules this work changed. Comparing a whole-tree run to a three-module
# baseline reports every unrelated pre-existing failure as a regression.
set +e
PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" -m unittest -v \
  tests.unit.test_check_contract_suite_baseline \
  tests.unit.test_config_write_transaction \
  tests.contract.test_registry_files \
  tests.unit.test_auprc_reference \
  tests.unit.test_stats \
  > "${OUT}/unittest.txt" 2>&1
set -e

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" scripts/check_contract_suite_baseline.py \
  --unittest-output "${OUT}/unittest.txt" \
  --require-module test_registry_files | tee "${OUT}/contract_baseline.txt"

tail -5 "${OUT}/unittest.txt"
echo "SUITE OUTPUT: ${OUT}"
