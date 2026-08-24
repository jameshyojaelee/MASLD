#!/usr/bin/env bash
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${OUTER_FOLD:=0}"
: "${LINEAGE:=hepatocyte}"
: "${SEED:=20260824}"
case "${OUTER_FOLD}" in 0|1|2|3|4) ;; *) printf 'Invalid OUTER_FOLD\n' >&2; exit 2 ;; esac
case "${LINEAGE}" in cholangiocyte|fibroblast|hepatocyte|macrophage|t_cell) ;;
  *) printf 'Invalid LINEAGE\n' >&2; exit 2 ;;
esac
case "${SEED}" in 20260824|20260825|20260826) ;;
  *) printf 'Invalid screening SEED\n' >&2; exit 2 ;;
esac
export PYTHONHASHSEED="${SEED}"
SPLIT_ID="donor${OUTER_FOLD}_genomic${OUTER_FOLD}"

if (( $# != 7 )); then
  printf 'usage: %s MODEL_ID PSEUDOBULK PSEUDOBULK_SHA INPUTS INPUTS_SHA BATCH_SIZE WRAPPER\n' "$0" >&2
  exit 2
fi
MODEL_ID=$1
PSEUDOBULK=$2
PSEUDOBULK_ARTIFACTS_SHA256=$3
INPUTS=$4
INPUTS_ARTIFACTS_SHA256=$5
BATCH_SIZE=$6
WRAPPER=$7
if [[ "${INPUTS_ARTIFACTS_SHA256}" == AUTO ]]; then
  if [[ ! -f "${INPUTS}/ARTIFACTS.json" ]]; then
    printf 'Dependency-complete input receipt is missing: %s\n' "${INPUTS}/ARTIFACTS.json" >&2
    exit 1
  fi
  read -r INPUTS_ARTIFACTS_SHA256 _ < <(sha256sum "${INPUTS}/ARTIFACTS.json")
fi
if [[ "${MODEL_ID}" != sequence_cnn_control && \
      "${MODEL_ID}" != sequence_transformer_control ]]; then
  printf 'Unrecognized sequence control: %s\n' "${MODEL_ID}" >&2
  exit 2
fi
if [[ ! "${BATCH_SIZE}" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Invalid batch size: %s\n' "${BATCH_SIZE}" >&2
  exit 2
fi
if [[ "${WRAPPER}" == *..* || -L "${WRAPPER}" || ! -f "${WRAPPER}" ]]; then
  printf 'Invalid production wrapper: %s\n' "${WRAPPER}" >&2
  exit 2
fi

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
FASTA_ROOT="${ROOT}/executions/sequence-fasta-21063808"
FASTA="${FASTA_ROOT}/GRCh38.p14.sequence_model.fa"
CONFIG="${ROOT}/config/sequence_control_architectures.json"
ARCHITECTURE="${ROOT}/scripts/sequence_control_architectures.py"
SMOKE="${ROOT}/executions/sequence-control-training-inference-probe-21064208"
UNIT="${ROOT}/executions/sequence-control-unit-21064228"
BATCH_PROBE="${ROOT}/executions/sequence-control-batch32-probe-21064243"
THROUGHPUT_PROBE="${ROOT}/executions/sequence-control-throughput-probe-21064254"
DETERMINISM="${ROOT}/executions/sequence-control-retrain-determinism-21064284"
TARGET="${ROOT}/executions/${MODEL_ID}-${LINEAGE}-${SPLIT_ID}-seed${SEED}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
PYTHONPATH_VALUE="${ROOT}:${CORE}/site-packages:${SOURCE}"
BIGWIG_RELATIVE="donor_test_fold_${OUTER_FOLD}/${LINEAGE}.tn5.bw"
BIGWIG="${PSEUDOBULK}/${BIGWIG_RELATIVE}"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed sequence-control campaign: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite sequence-control campaign: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/model" \
  "${STAGE}/predictions" \
  "${STAGE}/prepared" \
  "${STAGE}/validation"

printf '%s  %s\n' "${PSEUDOBULK_ARTIFACTS_SHA256}" \
  "${PSEUDOBULK}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${INPUTS_ARTIFACTS_SHA256}" \
  "${INPUTS}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 04bd59c825874f2d2d7dbf38b741633ccacad4735e95997a7477288e9c9656c4 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  "${FASTA_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' ed20dceeb315cffa2bc3de050d846eceefc552451fef741f0dd95b0149d10fce \
  "${FASTA}" | sha256sum --check --strict
printf '%s  %s\n' 9fc50e2fd96285e264c671f603dbcb4e53da44f0966608b0c4f274bc5d958b7b \
  "${SMOKE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 305a2e764d1e9c4d5acc0654c9eb7fb3c24572082f37517385e518dc5cbccd6e \
  "${UNIT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 84ede48ade629db5d9cf2afdd9cc4ea15a5de8c2563d901f4557506e0aeef288 \
  "${BATCH_PROBE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' aa1a0cd6f877ab486e641d7223734e8f0351671ef3d72f0c38689288cd0fd3e1 \
  "${THROUGHPUT_PROBE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 9dde66f5837bf388104070a967eed1de436cbd23bd14cd700332ddae9c3d1149 \
  "${DETERMINISM}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${PSEUDOBULK}" "${INPUTS}" "${SPLIT}" "${FASTA_ROOT}" \
  "${SMOKE}" "${UNIT}" "${BATCH_PROBE}" "${THROUGHPUT_PROBE}" \
  "${DETERMINISM}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("input frozen-tree verification passed")
PY
python - "${PSEUDOBULK}/training_pseudobulk_manifest.tsv" \
  "${PSEUDOBULK}" "${STAGE}/validation/pseudobulk_source.json" \
  "${OUTER_FOLD}" "${LINEAGE}" <<'PY'
from pathlib import Path
import csv
import json
from hashlib import sha256
import sys

manifest = Path(sys.argv[1]).resolve(strict=True)
root = Path(sys.argv[2]).resolve(strict=True)
outer_fold = int(sys.argv[4])
lineage = sys.argv[5]
with manifest.open(encoding="utf-8", newline="") as handle:
    rows = [
        row
        for row in csv.DictReader(handle, delimiter="\t")
        if int(row["donor_test_fold"]) == outer_fold and row["lineage_id"] == lineage
    ]
if len(rows) != 1:
    raise SystemExit("sequence-control pseudobulk row differs")
row = rows[0]
path = (root / row["bigwig_path"]).resolve(strict=True)
if root not in path.parents or path.is_symlink():
    raise SystemExit("pseudobulk bigWig path differs")
digest = sha256()
with path.open("rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
if path.stat().st_size != int(row["bigwig_size_bytes"]) or digest.hexdigest() != row["bigwig_sha256"]:
    raise SystemExit("pseudobulk bigWig bytes differ")
if (
    int(row["donor_valid_fold"]) != (outer_fold + 1) % 5
    or len({int(value) for value in row["donor_train_folds"].split(",")}) != 3
    or int(row["donors"]) < 1
    or int(row["tn5_insertions"]) != 2 * int(row["unique_fragments"])
):
    raise SystemExit("pseudobulk donor or Tn5 contract differs")
Path(sys.argv[3]).write_text(
    json.dumps(row, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY

module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env "PYTHONHASHSEED=${SEED}" \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python "${ROOT}/scripts/sequence_control_prepare_training.py" \
    --bigwig "${BIGWIG}" \
    --peaks "${INPUTS}/regions/chrombpnet.regions.bed" \
    --nonpeaks "${INPUTS}/nonpeaks/chrombpnet_negatives.bed" \
    --fold "${INPUTS}/regions/chrombpnet.fold.json" \
    --output "${STAGE}/prepared/data" \
    --seed "${SEED}" \
    --negative-sampling-ratio 0.1 \
    --outlier-threshold 0.9999 \
    --input-length 2114 \
    --output-length 1000 \
    --max-jitter 500 \
    > "${STAGE}/validation/preparation.stdout.txt" \
    2> "${STAGE}/validation/preparation.stderr.txt"

cd "${STAGE}/prepared/data"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env "PYTHONHASHSEED=${SEED}" \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python "${ROOT}/scripts/sequence_control_run_training.py" \
    --control-config "${CONFIG}" \
    --control-model-id "${MODEL_ID}" \
    --genome "${FASTA}" \
    --bigwig "${BIGWIG}" \
    --peaks sequence_control.filtered.peaks.bed \
    --nonpeaks sequence_control.filtered.nonpeaks.bed \
    --output_prefix "${STAGE}/model/${MODEL_ID}" \
    --chr_fold_path actual_fold.json \
    --epochs 50 \
    --early-stop 5 \
    --batch_size "${BATCH_SIZE}" \
    --learning-rate 0.001 \
    --params sequence_control_model_params.actual_fold.tsv \
    --seed "${SEED}" \
    --architecture_from_file "${ARCHITECTURE}" \
    > "${STAGE}/validation/training.stdout.txt" \
    2> "${STAGE}/validation/training.stderr.txt"

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env "PYTHONHASHSEED=${SEED}" \
  --env TF_DETERMINISTIC_OPS=1 \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
  "${IMAGE}" python "${ROOT}/scripts/sequence_control_predict_ccre.py" \
    --model-id "${MODEL_ID}" \
    --model "${STAGE}/model/${MODEL_ID}.h5" \
    --fasta "${FASTA}" \
    --windows "${SPLIT}/ccre_evaluation_windows.tsv" \
    --fold "${STAGE}/prepared/data/actual_fold.json" \
    --output "${STAGE}/predictions/fixed_ccre" \
    --lineage "${LINEAGE}" \
    --batch-size "${BATCH_SIZE}" \
    --seed "${SEED}" \
    --expected-windows-per-role 16000 \
    > "${STAGE}/validation/prediction.stdout.txt" \
    2> "${STAGE}/validation/prediction.stderr.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python - "${CONFIG}" "${INPUTS}" "${STAGE}" \
  "${MODEL_ID}" "${BATCH_SIZE}" "${LINEAGE}" "${SPLIT_ID}" "${SEED}" <<'PY'
from pathlib import Path
import csv
import json
import sys

import h5py
import numpy as np
import tensorflow as tf

from chrombpnet.training.utils.losses import multinomial_nll

config = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
inputs = Path(sys.argv[2]).resolve(strict=True)
root = Path(sys.argv[3]).resolve(strict=True)
model_id = sys.argv[4]
batch_size = int(sys.argv[5])
lineage = sys.argv[6]
split_id = sys.argv[7]
seed = int(sys.argv[8])
input_summary = json.loads((inputs / "validation" / "summary.json").read_text())
if (
    input_summary["split_id"] != split_id
    or input_summary["lineage_id"] != lineage
    or input_summary["test_outcomes_available_to_training"]
):
    raise SystemExit("production input summary differs")
preparation = json.loads((root / "prepared" / "data" / "summary.json").read_text())
if (
    preparation["hyperparameter_fit_outcome_roles"] != ["train"]
    or preparation["validation_outcomes_used_for_hyperparameter_fit"]
    or preparation["test_outcomes_used_for_hyperparameter_fit"]
    or preparation["bias_model_loaded"]
    or preparation["filtered_peak_role_census"]["valid"] != 16000
    or preparation["filtered_peak_role_census"]["test"] != 16000
):
    raise SystemExit("train-only preparation contract differs")
model = tf.keras.models.load_model(
    root / "model" / f"{model_id}.h5",
    custom_objects={"multinomial_nll": multinomial_nll, "tf": tf},
    compile=False,
)
budget = config["parameter_budget"]
if (
    model.name != model_id
    or model.input_shape != (None, 2114, 4)
    or [tuple(value) for value in model.output_shape] != [(None, 1000), (None, 1)]
    or not budget["minimum_trainable_parameters"]
    <= model.count_params()
    <= budget["maximum_trainable_parameters"]
):
    raise SystemExit("trained production model contract differs")
with (root / "model" / f"{model_id}.log").open(
    encoding="utf-8", newline=""
) as handle:
    history = list(csv.DictReader(handle))
if not 1 <= len(history) <= 50 or any(
    not np.isfinite(float(row["loss"])) or not np.isfinite(float(row["val_loss"]))
    for row in history
):
    raise SystemExit("training history differs")
receipt = json.loads(
    (root / "model" / f"{model_id}.training_receipt.json").read_text()
)
if receipt["test_outcomes_used"] or receipt["held_donor_atac_exposed"]:
    raise SystemExit("trained model receipt exposes held outcomes")
prediction = json.loads(
    (root / "predictions" / "fixed_ccre" / "summary.json").read_text()
)
contract = json.loads(
    (root / "predictions" / "fixed_ccre" / "prediction_contract.json").read_text()
)
if (
    prediction["role_counts"] != {"valid": 16000, "test": 16000}
    or prediction["observed_atac_input_exposed"]
    or prediction["benchmark_metrics_calculated"]
    or not contract["ready_for_independent_evaluator_expansion"]
):
    raise SystemExit("production prediction contract differs")
with h5py.File(
    root / "predictions" / "fixed_ccre" / "profile_probabilities.h5", "r"
) as handle:
    profiles = handle["profile_probability"]
    if profiles.shape != (32000, 1000):
        raise SystemExit("production profile geometry differs")
    for start in range(0, profiles.shape[0], 1024):
        values = profiles[start : start + 1024]
        if not np.isfinite(values).all() or not np.allclose(
            values.sum(axis=1), 1.0, rtol=0, atol=1e-6
        ):
            raise SystemExit("production profile values differ")
validation = {
    "schema_version": "masld-bench-sequence-control-production-v1",
    "status": "pass",
    "model_id": model_id,
    "dataset_id": "gse296875",
    "lineage_id": lineage,
    "split_id": split_id,
    "seed": seed,
    "batch_size": batch_size,
    "maximum_epochs": 50,
    "early_stopping_patience": 5,
    "epochs_completed": len(history),
    "parameter_count": model.count_params(),
    "input_length": 2114,
    "output_length": 1000,
    "fixed_ccre_windows": 32000,
    "reverse_complement_tta": True,
    "held_donor_atac_exposed": False,
    "test_outcomes_used": False,
    "benchmark_metrics_calculated": False,
    "ready_for_independent_evaluator_expansion": True,
}
(root / "validation" / "summary.json").write_text(
    json.dumps(validation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(validation, sort_keys=True))
PY
sha256sum \
  "${ROOT}/config/sequence_control_architectures.json" \
  "${ROOT}/scripts/sequence_control_architectures.py" \
  "${ROOT}/scripts/sequence_control_prepare_training.py" \
  "${ROOT}/scripts/sequence_control_run_training.py" \
  "${ROOT}/scripts/sequence_control_predict_ccre.py" \
  "${ROOT}/scripts/sequence_control_production_driver.sh" \
  "${WRAPPER}" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${PSEUDOBULK}" "${INPUTS}" "${SPLIT}" "${FASTA_ROOT}" \
  "${SMOKE}" "${UNIT}" "${BATCH_PROBE}" "${THROUGHPUT_PROBE}" \
  "${DETERMINISM}" \
  "${STAGE}/validation/input-verification.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

payload = {
    Path(value).name: verify_frozen_tree(Path(value)) for value in sys.argv[1:-1]
}
Path(sys.argv[-1]).write_text(
    json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${MODEL_ID}" \
  "${PSEUDOBULK_ARTIFACTS_SHA256}" "${INPUTS_ARTIFACTS_SHA256}" \
  "${BATCH_SIZE}" "${LINEAGE}" "${SPLIT_ID}" "${SEED}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "sequence_control_production",
            "model_id": sys.argv[2],
            "dataset_id": "gse296875",
            "lineage_id": sys.argv[6],
            "split_id": sys.argv[7],
            "pseudobulk_artifacts_sha256": sys.argv[3],
            "inputs_artifacts_sha256": sys.argv[4],
            "seed": int(sys.argv[8]),
            "batch_size": int(sys.argv[5]),
            "maximum_epochs": 50,
            "test_outcomes_used": False,
            "held_donor_atac_exposed": False,
            "status": "passed",
        },
    )
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Sequence-control production artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
