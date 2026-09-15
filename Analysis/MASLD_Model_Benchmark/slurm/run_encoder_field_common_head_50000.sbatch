#!/usr/bin/env bash
#SBATCH --job-name=model-encoder-field-7f
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Resource ask re-derived for this job, not copied from the scVI bundle:
#   memory   five frozen matrices 1.4G (1152+1280+1280+2048+128 float32 x 50k)
#            + fixture h5ad 0.6G + float64 normalized copy 0.9G
#            + HVG dense 50k x 2000 float32 0.4G + train dense 0.3G
#            + torch/CUDA host side ~2G  => ~6G peak, 48G requested (8x headroom).
#   time     210 shards x 6 head fits.  The frozen 30-shard Geneformer bundle
#            took 5:55 at 1152-d and TranscriptFormer 7:18 at 2048-d, so the
#            five frozen blocks are ~30 min; the two refit-per-fold blocks add
#            15 HVG+PCA passes and 15 latent loads; guards add ~15 min.
#            ~2 h expected, 12 h requested.
#   L40S     required: the frozen shards this job must reproduce were fit on
#            one L40S sm89 under CUBLAS_WORKSPACE_CONFIG=:4096:8.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONHASHSEED=20260824
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
REGISTRY=${ROOT}/config/encoder_field_common_head_50000_registry.json
TARGET=${ROOT}/executions/encoder-field-common-head-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  status=$?
  if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${status}"
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite the encoder-field bundle.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in \
  "${RUNTIME_SHA256} ${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE_SHA256} ${SOURCE}/ARTIFACTS.json" \
  "${SPLIT_SHA256} ${SPLIT}/ARTIFACTS.json" ; do
  printf '%s  %s\n' ${pair} | sha256sum --check --strict
done
export PYTHONPATH="${ROOT}:${ROOT}/src"

nvidia-smi -q > "${STAGE}/validation/nvidia-smi-q.txt"
"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"

# Guards first.  If a guard cannot be shown to fail on a deliberate tamper the
# job stops here and no score is produced.
"${PYTHON}" scripts/verify_encoder_field_common_head_50000.py \
  --source "${SOURCE}" \
  --split "${SPLIT}" \
  --registry "${REGISTRY}" \
  --output "${STAGE}/validation/guard_report.json" \
  > "${STAGE}/validation/guards.stdout.json" \
  2> "${STAGE}/validation/guards.stderr.txt"

sha256sum \
  scripts/fit_predict_encoder_field_common_head_50000.py \
  scripts/score_encoder_field_common_head_50000.py \
  scripts/verify_encoder_field_common_head_50000.py \
  scripts/fit_predict_common_cell_heads_study_50000.py \
  scripts/score_cell_baselines_study_50000.py \
  "${REGISTRY}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

# Shards, written one at a time.  A wall-clock timeout costs only the shard in
# flight; a resubmission skips every shard whose receipt already exists.
"${PYTHON}" scripts/fit_predict_encoder_field_common_head_50000.py \
  --source "${SOURCE}" \
  --split "${SPLIT}" \
  --registry "${REGISTRY}" \
  --output "${STAGE}" \
  > "${STAGE}/validation/shards.stdout.jsonl" \
  2> "${STAGE}/validation/shards.stderr.txt"

# Second reproduction control: sixty of these shards were already produced by
# the frozen common-head bundles.  They must come back identical.
"${PYTHON}" scripts/audit_encoder_field_shard_agreement.py \
  --shards "${STAGE}/shards" \
  --output "${STAGE}/validation/frozen_shard_agreement.json" \
  > "${STAGE}/validation/shard_agreement.stdout.json" \
  2> "${STAGE}/validation/shard_agreement.stderr.txt"

"${PYTHON}" - "${STAGE}" <<'PY'
import sys
from pathlib import Path

sys.path[:0] = [
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/src",
]
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
shards = sorted(path.name for path in (root / "shards").iterdir() if path.is_dir())
if len(shards) != 210:
    raise SystemExit(f"expected 210 shards, found {len(shards)}")
freeze_tree(
    root,
    {
        "artifact_class": "encoder_field_common_head_predictions",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "split_id": "resource_atlas_study_outer_5fold_v1",
        "shards": len(shards),
        "metrics_calculated": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
printf 'encoder field bundle: %s\n' "${TARGET}"
