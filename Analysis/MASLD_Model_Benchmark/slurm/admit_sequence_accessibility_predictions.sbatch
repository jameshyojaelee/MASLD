#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-prediction-admit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
MANIFEST="${ROOT}/config/evaluation/sequence_accessibility_d0g0_predictions.tsv"
cd "${ROOT}"
printf '%s  %s\n' 59ddf83b1f8206611234a8b28fa143dacbe900efadd2e539291f7eb0012c165d \
  config/evaluation/sequence_accessibility_d0g0_predictions.tsv \
  | sha256sum --check --strict
printf '%s  %s\n' ea2218ea75ad8d4f3c8097f40c0a37d19c925c9811d4342da0cee21e71e69302 \
  scripts/evaluate_sequence_accessibility.py | sha256sum --check --strict
printf '%s  %s\n' ec9806e1b613cdcd4baee0f73a1878f34c14fe018bbc9268233379f7cb6c28a7 \
  tests/unit/test_evaluate_sequence_accessibility.py \
  | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src" python - <<'PY'
from pathlib import Path
import csv

from masld_bench.artifacts import verify_frozen_tree

manifest = Path("config/evaluation/sequence_accessibility_d0g0_predictions.tsv")
with manifest.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
for row in rows:
    verify_frozen_tree(Path(row["root"]))
print(f"verified {len(rows)} frozen prediction roots")
PY
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2
PYTHONPATH="src:." python -m unittest \
  tests.unit.test_evaluate_sequence_accessibility -v
PYTHONPATH="src:." python - "${MANIFEST}" <<'PY'
from pathlib import Path
import sys

from scripts.evaluate_sequence_accessibility import (
    load_prediction,
    load_prediction_manifest,
)

rows = load_prediction_manifest(Path(sys.argv[1]))
for row in rows:
    prediction = load_prediction(row)
    print(
        row["model_id"],
        len(prediction["counts"]),
        prediction["profiles"].shape,
    )
PY
