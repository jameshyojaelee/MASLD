#!/usr/bin/env bash
#SBATCH --job-name=model-test-402
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
export PYTHONPATH="${root}/src:${root}"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "${root}"
module load python/3.11.5-GCCcore-13.2.0
python -m unittest -v tests.unit.test_gpu_bundle_dispatcher
