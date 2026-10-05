#!/bin/bash
#SBATCH --job-name=1596
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1596_h3_verify_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1596_h3_verify_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ $# -ne 1 ]]; then
    echo 'Require reviewed verifier source SHA256' >&2
    exit 2
fi
TASK_PRIMARY_SUMMARY_SHA256=$(/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python - <<'PY'
import hashlib
import json
from pathlib import Path

directory = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/codex_receiving_profile_evaluation_22026740')
raw = (directory / 'summary.json').read_bytes()
summary = json.loads(raw)
if ((directory / 'failure.json').exists()
        or summary.get('status') != 'complete_fixed_transport_evaluation'
        or summary.get('biological_n') != 17 or summary.get('regions') != 96460
        or summary.get('inputs_unchanged_after_evaluation') is not True
        or summary.get('model_fitted') is not False
        or summary.get('output_sha256', {}).get('executed_sources/codex_receiving_profile_evaluation.py')
        != '7ae9d401789880bd8ca40cbdbe75b54d7e73f43d1008048f6f390ee6308c9895'):
    raise SystemExit('Require completed exact frozen primary evaluation')
print(hashlib.sha256(raw).hexdigest())
PY
)
exec /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python codex_verify_receiving_profile_losses.py --expected-source-sha256 "$1" --expected-primary-summary-sha256 "$TASK_PRIMARY_SUMMARY_SHA256"
