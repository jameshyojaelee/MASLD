#!/usr/bin/env bash
#SBATCH --job-name=model-test-304
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
runtime="${root}/executions/environments/scprint2-native-torch2.8.0-v2-21073597/source/.venv"
export PYTHONPATH="${root}/src"
export LD_LIBRARY_PATH="/gpfs/commons/home/jameslee/micromamba/envs/scprint2x/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
cd "${root}"
bash -n slurm/model_work_302.sbatch
bash -n slurm/evaluate_scprint2_common_lane.sbatch
"${runtime}/bin/python" -m py_compile scripts/scprint2_extract_no_neighbor_embeddings.py
/gpfs/commons/home/jameslee/micromamba/envs/gca-test/bin/python -m py_compile \
  scripts/evaluate_transcriptformer_common_lane.py
/gpfs/commons/home/jameslee/micromamba/envs/gca-test/bin/python -m pytest -q \
  tests/unit/test_scprint2_embedding_contract.py
