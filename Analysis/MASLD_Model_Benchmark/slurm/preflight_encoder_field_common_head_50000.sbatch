#!/usr/bin/env bash
#SBATCH --job-name=model-encoder-field-7f-preflight
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# CPU-only preflight.  Exercises every code path that does not need the GPU:
# block construction, both refit-per-fold feature builders, and the four CPU
# guards including the scVI latent-placement reproduction.  Cheap insurance
# against burning GPU hours on a typo.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONHASHSEED=20260824
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SPLIT=${ROOT}/executions/model-data-052-21077595/split
REGISTRY=${ROOT}/config/encoder_field_common_head_50000_registry.json
TARGET=${ROOT}/executions/encoder-field-preflight-${SLURM_JOB_ID}

install -d -m 0750 "${TARGET}"
cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"

"${PYTHON}" - "${SOURCE}" "${SPLIT}" "${REGISTRY}" "${TARGET}" <<'PY'
import json
import sys
from pathlib import Path

import numpy as np

from scripts.fit_predict_encoder_field_common_head_50000 import build_blocks, load_split

source, split, registry_path, target = (Path(value) for value in sys.argv[1:5])
registry = json.loads(registry_path.read_text(encoding="utf-8"))
context = load_split(source, split)
records = []
for block in build_blocks(registry=registry, context=context, fixture=source):
    for held_outer in (0, 4):
        features = block.features(screen_seed=1103, held_outer=held_outer)
        if features.shape != (50_000, block.width) or not np.all(np.isfinite(features)):
            raise SystemExit(f"{block.block_id}: fold {held_outer} features differ")
        train = context["outer"] != held_outer
        records.append(
            {
                "block_id": block.block_id,
                "family": block.family,
                "outer_fold": held_outer,
                "width": block.width,
                "training_rows": int(train.sum()),
                "held_out_rows": int((~train).sum()),
                "feature_abs_mean": float(np.abs(features).mean()),
                "training_feature_variance_mean": float(features[train].var(axis=0).mean()),
            }
        )
        print(json.dumps(records[-1], sort_keys=True), flush=True)
(target / "feature_preflight.json").write_text(
    json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY

"${PYTHON}" scripts/verify_encoder_field_common_head_50000.py \
  --source "${SOURCE}" \
  --split "${SPLIT}" \
  --registry "${REGISTRY}" \
  --output "${TARGET}/guard_report_cpu_only.json" \
  --skip-gpu-checks

printf 'encoder field preflight: %s\n' "${TARGET}"
