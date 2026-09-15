#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-606
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
export PYTHONPATH=${ROOT}:${ROOT}/src
cd "${ROOT}"

"${PYTHON}" -m unittest -v \
  tests.unit.test_gse281364_dna_language_seeded_campaign

"${PYTHON}" scripts/check_gse281364_shared_control_copy_preflight.py \
  --root "${ROOT}" \
  --config "${ROOT}/config/gse281364_dna_language_seeded_head_campaign.json"

"${PYTHON}" -m json.tool \
  config/gse281364_dna_language_seeded_attempt_dispositions.json \
  > /dev/null
