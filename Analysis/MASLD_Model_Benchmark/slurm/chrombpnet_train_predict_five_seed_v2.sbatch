#!/usr/bin/env bash
#SBATCH --job-name=model-training-702
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${OUTER_FOLD:=0}"
: "${LINEAGE:=hepatocyte}"
: "${SEED:=20260824}"
: "${INPUT_JOB_ID:=21064244}"
: "${INPUT_ARTIFACTS_SHA256:=fc26ec558ea51f74212bea120dab59053e5f2856494407bb0a717e9eaecc0208}"
case "${OUTER_FOLD}" in 0|1|2|3|4) ;; *) printf 'Invalid OUTER_FOLD\n' >&2; exit 2 ;; esac
case "${LINEAGE}" in cholangiocyte|fibroblast|hepatocyte|macrophage|t_cell) ;;
  *) printf 'Invalid LINEAGE\n' >&2; exit 2 ;;
esac
case "${SEED}" in 20260824|20260825|20260826|20260827|20260828) ;;
  *) printf 'Invalid screening SEED\n' >&2; exit 2 ;;
esac
export PYTHONHASHSEED="${SEED}"
SPLIT_ID="donor${OUTER_FOLD}_genomic${OUTER_FOLD}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
BIAS_ROOT="${ROOT}/executions/chrombpnet-bias-models-zenodo7443683-21063957"
BIAS_MODEL="${BIAS_ROOT}/models/bias_models/ATAC/scATAC_dermal_fibroblast.h5"
PSEUDOBULK="${ROOT}/executions/gse296875-training-fold-pseudobulk-21063993"
PSEUDOBULK_ARTIFACTS_SHA256=3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b
TRAIN_BIGWIG="${PSEUDOBULK}/donor_test_fold_${OUTER_FOLD}/${LINEAGE}.tn5.bw"
INPUT_ROOT="${ROOT}/executions/chrombpnet-${LINEAGE}-${SPLIT_ID}-inputs-${INPUT_JOB_ID}"
if [[ "${LINEAGE}:${OUTER_FOLD}:${INPUT_JOB_ID}" == hepatocyte:0:21064244 ]]; then
  INPUT_ROOT="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-inputs-21064244"
fi
if [[ "${INPUT_ARTIFACTS_SHA256}" == AUTO ]]; then
  if [[ ! -f "${INPUT_ROOT}/ARTIFACTS.json" ]]; then
    printf 'Dependency-complete input receipt is missing: %s\n' "${INPUT_ROOT}/ARTIFACTS.json" >&2
    exit 1
  fi
  read -r INPUT_ARTIFACTS_SHA256 _ < <(sha256sum "${INPUT_ROOT}/ARTIFACTS.json")
fi
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
WINDOWS="${SPLIT}/ccre_evaluation_windows.tsv"
FASTA_ROOT="${ROOT}/executions/sequence-fasta-21063808"
FASTA="${FASTA_ROOT}/GRCh38.p14.sequence_model.fa"
TARGET="${ROOT}/executions/chrombpnet-${LINEAGE}-${SPLIT_ID}-seed${SEED}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
PYTHONPATH_VALUE="${ROOT}:${CORE}/site-packages:${SOURCE}"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed ChromBPNet training stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite ChromBPNet training artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
if [[ "${PSEUDOBULK_ARTIFACTS_SHA256}" == TO_BE_FILLED* \
  || "${INPUT_ROOT}" == *TO_BE_FILLED* \
  || "${INPUT_ARTIFACTS_SHA256}" == TO_BE_FILLED* ]]; then
  printf 'Production input bindings are unresolved.\n' >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/prepared" \
  "${STAGE}/model" \
  "${STAGE}/predictions" \
  "${STAGE}/validation"

printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 04bd59c825874f2d2d7dbf38b741633ccacad4735e95997a7477288e9c9656c4 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 8ecc40f452328df965be823d14f9ce76e2777407dba62cc7826685be635bb6e5 \
  "${BIAS_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 3fcadc521716156405650b1f439bce81c074193c226b4883471b328f1c7b27be \
  "${BIAS_MODEL}" | sha256sum --check --strict
printf '%s  %s\n' "${PSEUDOBULK_ARTIFACTS_SHA256}" \
  "${PSEUDOBULK}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${INPUT_ARTIFACTS_SHA256}" \
  "${INPUT_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 3dbefbb961c040236431773fe3e2e32527622e250e9fc86133ed117dc4420c4c \
  "${WINDOWS}" | sha256sum --check --strict
printf '%s  %s\n' 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  "${FASTA_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' ed20dceeb315cffa2bc3de050d846eceefc552451fef741f0dd95b0149d10fce \
  "${FASTA}" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${PSEUDOBULK}" "${INPUT_ROOT}" "${SPLIT}" "${FASTA_ROOT}" "${BIAS_ROOT}" \
  "${STAGE}/validation/input-verification.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

roots = [Path(value).resolve(strict=True) for value in sys.argv[1:-1]]
if any("outcome" in part.lower() for root in roots for part in root.parts):
    raise SystemExit("evaluator outcome artifact entered the model job")
payload = {root.name: verify_frozen_tree(root) for root in roots}
Path(sys.argv[-1]).write_text(
    json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print("model-input frozen-tree and outcome-firewall verification passed")
PY

module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
cd "${ROOT}"

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env "PYTHONHASHSEED=${SEED}" \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python scripts/prepare_chrombpnet_train_only_hyperparams.py \
    --genome "${FASTA}" \
    --bigwig "${TRAIN_BIGWIG}" \
    --peaks "${INPUT_ROOT}/regions/chrombpnet.regions.bed" \
    --nonpeaks "${INPUT_ROOT}/nonpeaks/chrombpnet_negatives.bed" \
    --actual-fold "${INPUT_ROOT}/regions/chrombpnet.fold.json" \
    --bias-model "${BIAS_MODEL}" \
    --output "${STAGE}/prepared/native" \
    --seed "${SEED}" \
    --negative-sampling-ratio 0.1 \
    --outlier-threshold 0.9999 \
    --input-length 2114 \
    --output-length 1000 \
    --max-jitter 500 \
    --filters 512 \
    --dilation-layers 8 \
    > "${STAGE}/validation/hyperparameter-preparation.stdout.txt" \
    2> "${STAGE}/validation/hyperparameter-preparation.stderr.txt"

cd "${STAGE}/prepared/native"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env "PYTHONHASHSEED=${SEED}" \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python "${ROOT}/scripts/run_chrombpnet_training.py" \
    --genome "${FASTA}" \
    --bigwig "${TRAIN_BIGWIG}" \
    --peaks train_only_filtered.peaks.bed \
    --nonpeaks train_only_filtered.nonpeaks.bed \
    --output_prefix "${STAGE}/model/chrombpnet" \
    --chr_fold_path actual_fold.json \
    --epochs 50 \
    --early-stop 5 \
    --batch_size 64 \
    --learning-rate 0.001 \
    --params chrombpnet_model_params.actual_fold.tsv \
    --seed "${SEED}" \
    --architecture_from_file \
      "${SOURCE}/chrombpnet/training/models/chrombpnet_with_bias_model.py" \
    > "${STAGE}/validation/training.stdout.txt" \
    2> "${STAGE}/validation/training.stderr.txt"

cd "${ROOT}"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env "PYTHONHASHSEED=${SEED}" \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python scripts/predict_chrombpnet_ccre_windows.py \
    --model "${STAGE}/model/chrombpnet.h5" \
    --fasta "${FASTA}" \
    --windows "${WINDOWS}" \
    --fold "${INPUT_ROOT}/regions/chrombpnet.fold.json" \
    --output "${STAGE}/predictions/full_model" \
    --batch-size 512 \
    --seed "${SEED}" \
    --count-transform log1p_absolute \
    > "${STAGE}/validation/full-model-prediction.stdout.txt" \
    2> "${STAGE}/validation/full-model-prediction.stderr.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env "PYTHONHASHSEED=${SEED}" \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python scripts/predict_chrombpnet_ccre_windows.py \
    --model "${STAGE}/model/chrombpnet_nobias.h5" \
    --fasta "${FASTA}" \
    --windows "${WINDOWS}" \
    --fold "${INPUT_ROOT}/regions/chrombpnet.fold.json" \
    --output "${STAGE}/predictions/nobias_model" \
    --batch-size 512 \
    --seed "${SEED}" \
    --count-transform log_component \
    > "${STAGE}/validation/nobias-model-prediction.stdout.txt" \
    2> "${STAGE}/validation/nobias-model-prediction.stderr.txt"

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  "${IMAGE}" python - "${STAGE}" "${INPUT_ROOT}" \
    "${LINEAGE}" "${SPLIT_ID}" "${SEED}" <<'PY'
from pathlib import Path
import csv
from hashlib import sha256
import json
import sys

import h5py
import numpy as np
import tensorflow as tf

from chrombpnet.training.utils.losses import multinomial_nll

root = Path(sys.argv[1]).resolve(strict=True)
inputs = Path(sys.argv[2]).resolve(strict=True)
lineage = sys.argv[3]
split_id = sys.argv[4]
seed = int(sys.argv[5])
prepared = json.loads((root / "prepared/native/summary.json").read_text())
input_summary = json.loads((inputs / "validation/summary.json").read_text())
if (
    input_summary["lineage_id"] != lineage
    or input_summary["split_id"] != split_id
    or input_summary["test_outcomes_available_to_training"]
    or prepared["hyperparameter_fit_outcome_roles"] != ["train"]
    or prepared["validation_outcomes_used_for_hyperparameter_fit"]
    or prepared["test_outcomes_used_for_hyperparameter_fit"]
    or prepared["filters"] != 512
    or prepared["dilation_layers"] != 8
):
    raise SystemExit("train-only preparation contract differs")
custom_objects = {"multinomial_nll": multinomial_nll, "tf": tf}
model_paths = {
    "full_model": root / "model/chrombpnet.h5",
    "nobias_model": root / "model/chrombpnet_nobias.h5",
}
model_sha256 = {}
for model_id, path in model_paths.items():
    if not path.is_file():
        raise SystemExit(f"missing trained checkpoint: {model_id}")
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    model_sha256[model_id] = digest.hexdigest()
    model = tf.keras.models.load_model(path, custom_objects=custom_objects, compile=False)
    if model.input_shape != (None, 2114, 4) or [tuple(value) for value in model.output_shape] != [
        (None, 1000),
        (None, 1),
    ]:
        raise SystemExit(f"trained model topology differs: {model_id}")
with (root / "model/chrombpnet.log").open(encoding="utf-8", newline="") as handle:
    history = list(csv.DictReader(handle))
if not 1 <= len(history) <= 50 or any(
    not np.isfinite(float(row["loss"])) or not np.isfinite(float(row["val_loss"]))
    for row in history
):
    raise SystemExit("training history differs")
arguments = json.loads((root / "model/chrombpnet.args.json").read_text())
if (
    arguments["seed"] != seed
    or arguments["epochs"] != 50
    or arguments["early_stop"] != 5
    or arguments["batch_size"] != 64
):
    raise SystemExit("captured training arguments differ")
for model_id, transform in {
    "full_model": "log1p_absolute",
    "nobias_model": "log_component",
}.items():
    summary = json.loads((root / f"predictions/{model_id}/summary.json").read_text())
    if (
        summary["role_counts"] != {"valid": 16000, "test": 16000}
        or summary["regional_mass_transform"] != transform
        or summary["observed_atac_input_exposed"]
        or summary["benchmark_metrics_calculated"]
    ):
        raise SystemExit(f"prediction firewall differs: {model_id}")
    with h5py.File(root / f"predictions/{model_id}/profile_probabilities.h5", "r") as handle:
        values = handle["profile_probability"]
        if values.shape != (32000, 1000):
            raise SystemExit(f"prediction tensor shape differs: {model_id}")
validation = {
    "schema_version": "masld-bench-chrombpnet-full-depth-training-v1",
    "status": "pass",
    "dataset_id": "gse296875",
    "lineage_id": lineage,
    "split_id": split_id,
    "seed": seed,
    "filters": 512,
    "dilation_layers": 8,
    "epochs_requested": 50,
    "epochs_completed": len(history),
    "validation_outcomes_used_for_early_stopping": True,
    "test_outcomes_used": False,
    "observed_held_donor_atac_used_for_inference": False,
    "model_sha256": model_sha256,
    "primary_biological_score": "nobias_model",
    "assay_qc_score": "full_model",
    "source_bias_role": "restricted_smoke_comparator",
    "champion_eligible": False,
    "benchmark_metrics_calculated": False,
}
(root / "validation/summary.json").write_text(
    json.dumps(validation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(validation, sort_keys=True))
PY

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/prepare_chrombpnet_train_only_hyperparams.py" \
  "${ROOT}/scripts/run_chrombpnet_training.py" \
  "${ROOT}/scripts/predict_chrombpnet_ccre_windows.py" \
  "${ROOT}/slurm/chrombpnet_train_predict_five_seed_v2.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" \
  "${PSEUDOBULK_ARTIFACTS_SHA256}" "${INPUT_ARTIFACTS_SHA256}" \
  "${LINEAGE}" "${SPLIT_ID}" "${SEED}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "chrombpnet_full_depth_training_and_predictions",
            "model_id": "chrombpnet",
            "source_revision": "eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c",
            "dataset_id": "gse296875",
            "lineage_id": sys.argv[4],
            "split_id": sys.argv[5],
            "seed": int(sys.argv[6]),
            "pseudobulk_artifacts_sha256": sys.argv[2],
            "training_input_artifacts_sha256": sys.argv[3],
            "filters": 512,
            "dilation_layers": 8,
            "epochs_requested": 50,
            "early_stop": 5,
            "batch_size": 64,
            "primary_biological_score": "nobias_model",
            "assay_qc_score": "full_model",
            "source_bias_sha256": "3fcadc521716156405650b1f439bce81c074193c226b4883471b328f1c7b27be",
            "source_bias_role": "restricted_smoke_comparator",
            "champion_eligible": False,
            "test_outcomes_used": False,
            "benchmark_metrics_calculated": False,
            "status": "passed",
        },
    )
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'ChromBPNet training/prediction campaign: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
