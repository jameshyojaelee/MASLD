#!/usr/bin/env bash
#SBATCH --job-name=masld-seq-training-mean
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
PSEUDOBULK="${ROOT}/executions/gse296875-training-fold-pseudobulk-21063993"
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
TARGET="${ROOT}/executions/sequence-training-mean-hepatocyte-donor0-genomic0-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite sequence training-mean baseline: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b \
  "${PSEUDOBULK}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${PSEUDOBULK}" "${SPLIT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("input frozen-tree verification passed")
PY
module purge
module load singularity/4.2.2
cd "${ROOT}"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python -m unittest -v \
    tests.unit.test_build_sequence_training_mean_baseline \
    > "${STAGE}/validation/tests.stdout.txt" \
    2> "${STAGE}/validation/tests.stderr.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python scripts/build_sequence_training_mean_baseline.py \
    --pseudobulk "${PSEUDOBULK}" \
    --pseudobulk-artifacts-sha256 \
      3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b \
    --split-contract "${SPLIT}" \
    --split-artifacts-sha256 \
      00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
    --split-id donor0_genomic0 \
    --lineage hepatocyte \
    --output "${STAGE}/baseline" \
    > "${STAGE}/validation/build.stdout.txt" \
    2> "${STAGE}/validation/build.stderr.txt"
singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${ROOT}:${CORE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

import h5py

root = Path(sys.argv[1]).resolve(strict=True)
summary = json.loads((root / "baseline" / "summary.json").read_text())
if (
    summary["roles"] != {"test": 16000, "valid": 16000}
    or summary["training_donor_folds"] != "2,3,4"
    or summary["training_donors"] != 19
    or summary["held_donor_outcomes_exposed"]
    or summary["benchmark_metrics_calculated"]
):
    raise SystemExit("training-mean summary differs")
with h5py.File(root / "baseline" / "training_mean_baseline.h5", "r") as handle:
    for role in ("valid", "test"):
        if handle[role]["pseudobulk_counts"].shape != (16000, 1000):
            raise SystemExit(f"{role} baseline profile shape differs")
validation = {
    "schema_version": "masld-bench-sequence-training-mean-probe-v1",
    "status": "pass",
    "training_donor_folds": "2,3,4",
    "training_donors": 19,
    "windows": 32000,
    "held_donor_outcomes_exposed": False,
}
(root / "validation" / "probe.json").write_text(
    json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(validation, sort_keys=True))
PY
sha256sum \
  scripts/build_sequence_training_mean_baseline.py \
  tests/unit/test_build_sequence_training_mean_baseline.py \
  slurm/build_sequence_training_mean_baseline.sbatch \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "sequence_training_mean_baseline",
            "model_id": "training_pseudobulk_mean",
            "dataset_id": "gse296875",
            "split_id": "donor0_genomic0",
            "lineage_id": "hepatocyte",
            "biological_unit": "donor",
            "training_donor_folds": "2,3,4",
            "training_donors": 19,
            "held_donor_outcomes_exposed": False,
            "benchmark_metrics_calculated": False,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Sequence training-mean baseline: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
