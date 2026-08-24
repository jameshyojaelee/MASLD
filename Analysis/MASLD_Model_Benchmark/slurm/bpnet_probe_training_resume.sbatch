#!/usr/bin/env bash
#SBATCH --job-name=masld-bpnet-resume-probe
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONHASHSEED=20260824

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/bpnet-v2.0.0-21064186"
SOURCE="${ACQUISITION}/source/bpnet-f4593eedca51741f25b8c5cc0c0d647faf8c8a3c"
TRAINING_PROBE="${ROOT}/executions/bpnet-training-inference-probe-21064291"
FIXTURE="${TRAINING_PROBE}/fixture"
TARGET="${ROOT}/executions/bpnet-training-resume-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${CORE}/site-packages:${SOURCE}:${ROOT}/scripts"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite BPNet resume probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 8062574bfaa9067b0aca82ffa3eeae87188140ee5f5842399171ad4535cf88f6 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' ed43000954718d58870b55fcfed6ac9ae1fb0df1a853d777ac83c99cdb0ab051 \
  "${TRAINING_PROBE}/ARTIFACTS.json" | sha256sum --check --strict

module load singularity/4.2.2
run_training() {
  local output="$1"
  shift
  singularity exec --nv --cleanenv \
    --bind "${ROOT}:${ROOT}:rw" \
    --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
    --env PYTHONNOUSERSITE=1 \
    --env PYTHONDONTWRITEBYTECODE=1 \
    --env PYTHONHASHSEED=20260824 \
    --env TF_DETERMINISTIC_OPS=1 \
    --env TF_FORCE_GPU_ALLOW_GROWTH=true \
    --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
    "${IMAGE}" python "${ROOT}/scripts/bpnet_run_training.py" \
      --genome "${FIXTURE}/fixture.fa" \
      --chrom-sizes "${FIXTURE}/fixture.chrom.sizes" \
      --bigwig "${FIXTURE}/fixture.tn5.bw" \
      --regions "${FIXTURE}/regions.narrowPeak" \
      --nonpeaks "${FIXTURE}/nonpeaks.bed" \
      --fold "${FIXTURE}/fold.json" \
      --output "${output}" \
      --epochs 3 \
      --batch-size 4 \
      --threads 1 \
      --max-jitter 8 \
      --negative-ratio 0.5 \
      --early-stop 5 \
      --seed 20260824 \
      "$@"
}

run_training "${STAGE}/interrupted" --attempt-epoch-limit 1 \
  > "${STAGE}/validation/interrupted.stdout.txt" \
  2> "${STAGE}/validation/interrupted.stderr.txt"
run_training "${STAGE}/resumed" --resume-from "${STAGE}/interrupted" \
  > "${STAGE}/validation/resumed.stdout.txt" \
  2> "${STAGE}/validation/resumed.stderr.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  "${IMAGE}" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
interrupted = json.loads((root / "interrupted" / "training_receipt.json").read_text())
resumed = json.loads((root / "resumed" / "training_receipt.json").read_text())
if interrupted["completed_epochs"] != 1 or interrupted["campaign_complete"]:
    raise SystemExit("interrupted BPNet receipt differs")
if resumed["completed_epochs"] != 3 or not resumed["campaign_complete"]:
    raise SystemExit("resumed BPNet receipt differs")
if resumed["training_signature"] != interrupted["training_signature"]:
    raise SystemExit("resume training signature differs")
validation = resumed["resume_validation"]
if (
    validation is None
    or validation["restored_completed_epochs"] != 1
    or validation["restoration_prediction_max_abs"] > 1e-6
    or validation["restored_optimizer_iterations"] < 1
):
    raise SystemExit("BPNet restored-state validation differs")
interrupted_history = json.loads((root / "interrupted" / "history.json").read_text())
resumed_history = json.loads((root / "resumed" / "history.json").read_text())
if resumed_history[0] != interrupted_history[0] or len(resumed_history) != 3:
    raise SystemExit("BPNet resumed history did not preserve completed epoch")
(root / "validation" / "probe.json").write_text(
    json.dumps({
        "schema_version": "masld-bench-bpnet-training-resume-probe-v1",
        "status": "pass",
        "seed": 20260824,
        "interrupted_epochs": 1,
        "resumed_epochs": 3,
        "resume_prediction_max_abs": validation["restoration_prediction_max_abs"],
        "resume_atol": 1e-6,
        "evaluator_outcomes_exposed": False,
    }, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
print(validation["restoration_prediction_max_abs"])
PY

sha256sum \
  "${ROOT}/scripts/bpnet_contract.py" \
  "${ROOT}/scripts/bpnet_run_training.py" \
  "${ROOT}/slurm/bpnet_probe_training_resume.sbatch" \
  > "${STAGE}/validation/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(freeze_tree(root, {
    "artifact_class": "bpnet_training_resume_probe",
    "model_id": "bpnet",
    "seed": 20260824,
    "resume_prediction_atol": 1e-6,
    "evaluator_outcomes_exposed": False,
    "status": "passed",
}))
PY
mv "${STAGE}" "${TARGET}"
printf 'BPNet training resume probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
