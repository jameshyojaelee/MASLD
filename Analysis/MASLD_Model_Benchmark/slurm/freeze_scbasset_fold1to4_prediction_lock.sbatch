#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-lock-032
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-lock-032-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-cpu-lock-032-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BUILDER=${ROOT}/scripts/freeze_scbasset_crossfold_prediction_lock.py
SOURCE_LOCK=${ROOT}/executions/scbasset-crossfold-evaluator-source-21070472
CHAIN=${ROOT}/executions/scbasset-fold1to4-chain-21069969
OUTPUT=${ROOT}/executions/scbasset-fold1to4-prediction-lock-${SLURM_JOB_ID}
declare -a INPUTS=(
  "${ROOT}/executions/scbasset-donor1_genomic1-inputs-21067171"
  "${ROOT}/executions/scbasset-donor2_genomic2-inputs-21067172"
  "${ROOT}/executions/scbasset-donor3_genomic3-inputs-21067173"
  "${ROOT}/executions/scbasset-donor4_genomic4-inputs-21067174"
)

module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
printf '%s  %s\n' \
  be05e2b17d22de121880ca91b3b9750b0800f049731c05d1e243d7ad0441e9ad "${SOURCE_LOCK}/ARTIFACTS.json" \
  9c113b082efcb9108b36d93b644687f4877a661097faa369b4c436023551ec48 "${BUILDER}" \
  | sha256sum --check --strict
PYTHONPATH=${ROOT}/src python - "${SOURCE_LOCK}" "${CHAIN}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
arguments=()
for input_root in "${INPUTS[@]}"; do
  arguments+=(--input-root "${input_root}")
done
PYTHONPATH=${ROOT}/src python "${BUILDER}" \
  --chain-root "${CHAIN}" \
  "${arguments[@]}" \
  --output "${OUTPUT}"
PYTHONPATH=${ROOT}/src python - "${OUTPUT}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root)
lock = json.loads((root / "prediction_lock.json").read_text(encoding="utf-8"))
if (
    manifest["metadata"]["outcomes_read"] is not False
    or lock["folds"] != [1, 2, 3, 4]
    or lock["seed"] != 20260824
    or lock["outcomes_read"] is not False
    or lock["outcome_paths_recorded"] is not False
    or lock["champion_claim_allowed"] is not False
):
    raise SystemExit("scBasset prediction lock contract differs")
print(sha256_file(root / "ARTIFACTS.json"))
PY
