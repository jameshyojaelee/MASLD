#!/usr/bin/env bash
#SBATCH --job-name=model-test-504
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
python_bin=${root}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
export PYTHONPATH="${root}:${root}/src"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
cd "${root}"
"${python_bin}" -m unittest -v tests.unit.test_score_cell_baselines_study_50000
