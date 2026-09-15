#!/usr/bin/env bash
#SBATCH --job-name=model-check-347
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "${ROOT}"; module purge; module load h5py/3.9.0-foss-2023a
# Verification only: no staging, no freeze, no config write.
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_mpralegnet_deliberate_non_changes \
  tests.unit.test_mpralegnet_terminal_disposition
