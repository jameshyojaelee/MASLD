#!/usr/bin/env bash
#SBATCH --job-name=masld-sequence-control-probe
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
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
ACQUISITION="${ROOT}/executions/upstream-acquisitions/chrombpnet-v1.0.1-21063022"
SOURCE="${ACQUISITION}/source/chrombpnet-eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c"
INPUT="${ROOT}/executions/chrombpnet-train-only-hyperparams-probe-21064020"
UNIT="${ROOT}/executions/sequence-control-unit-21064202"
CONFIG="${ROOT}/config/sequence_control_architectures.json"
ARCHITECTURE="${ROOT}/scripts/sequence_control_architectures.py"
TARGET="${ROOT}/executions/sequence-control-training-inference-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
PYTHONPATH_VALUE="${ROOT}:${CORE}/site-packages:${SOURCE}"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite sequence-control probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/fixture" \
  "${STAGE}/models" \
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
printf '%s  %s\n' 5713097b91248632eb8486cd7e6c59d4c1f53832f69f5febca1631d113b480d2 \
  "${INPUT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 7b6723d720bf66bf00db078b0eae9bd82fb1761fd6901dc9cf70ba040d88ff5d \
  "${UNIT}/ARTIFACTS.json" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${INPUT}" "${UNIT}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("input frozen-tree verification passed")
PY
python - "${INPUT}/prepared/actual_fold.json" \
  "${STAGE}/fixture/mini_ccre_windows.tsv" <<'PY'
from pathlib import Path
import csv
from hashlib import sha256
import json
import sys

fold = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
fields = (
    "contig", "output_start", "output_end", "window_id", "genomic_fold",
    "window_class", "ccre_class", "ccre_id", "ccre_start", "ccre_end",
    "input_start", "input_end", "selection_hash",
)
rows = []
for role, genomic_fold in (("valid", 1), ("test", 0)):
    for index, contig in enumerate(fold[role][:4]):
        center = 4000 + 1000 * index
        output_start = center - 500
        input_start = center - 1057
        identity = f"{role}:{contig}:{center}"
        window_id = sha256(identity.encode("utf-8")).hexdigest()
        selection_hash = sha256((identity + ":selection").encode("utf-8")).hexdigest()
        rows.append(
            {
                "contig": contig,
                "output_start": output_start,
                "output_end": output_start + 1000,
                "window_id": window_id,
                "genomic_fold": genomic_fold,
                "window_class": "ccre",
                "ccre_class": "PLS",
                "ccre_id": f"probe_{role}_{index}",
                "ccre_start": center - 100,
                "ccre_end": center + 100,
                "input_start": input_start,
                "input_end": input_start + 2114,
                "selection_hash": selection_hash,
            }
        )
with Path(sys.argv[2]).open("x", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
PY

module purge
module load singularity/4.2.2
nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv,noheader \
  > "${STAGE}/validation/nvidia-smi.csv"
singularity --version > "${STAGE}/validation/singularity.version.txt"
cd "${INPUT}/prepared"
for model_id in sequence_cnn_control sequence_transformer_control; do
  singularity exec --nv --cleanenv \
    --bind "${ROOT}:${ROOT}:rw" \
    --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
    --env PYTHONNOUSERSITE=1 \
    --env PYTHONDONTWRITEBYTECODE=1 \
    --env PYTHONHASHSEED=20260824 \
    --env TF_DETERMINISTIC_OPS=1 \
    --env TF_FORCE_GPU_ALLOW_GROWTH=true \
    --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
    "${IMAGE}" python "${ROOT}/scripts/sequence_control_run_training.py" \
      --control-config "${CONFIG}" \
      --control-model-id "${model_id}" \
      --genome "${INPUT}/fixture/fixture.fa" \
      --bigwig "${INPUT}/fixture/fixture.tn5.bw" \
      --peaks train_only_filtered.peaks.bed \
      --nonpeaks train_only_filtered.nonpeaks.bed \
      --output_prefix "${STAGE}/models/${model_id}" \
      --chr_fold_path actual_fold.json \
      --epochs 1 \
      --early-stop 1 \
      --batch_size 4 \
      --learning-rate 0.001 \
      --params chrombpnet_model_params.actual_fold.tsv \
      --seed 20260824 \
      --architecture_from_file "${ARCHITECTURE}" \
      > "${STAGE}/validation/${model_id}.training.stdout.txt" \
      2> "${STAGE}/validation/${model_id}.training.stderr.txt"

  singularity exec --nv --cleanenv \
    --bind "${ROOT}:${ROOT}:rw" \
    --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
    --env PYTHONNOUSERSITE=1 \
    --env PYTHONDONTWRITEBYTECODE=1 \
    --env PYTHONHASHSEED=20260824 \
    --env TF_DETERMINISTIC_OPS=1 \
    --env TF_FORCE_GPU_ALLOW_GROWTH=true \
    --env CUBLAS_WORKSPACE_CONFIG=:4096:8 \
    "${IMAGE}" python "${ROOT}/scripts/sequence_control_predict_ccre.py" \
      --model-id "${model_id}" \
      --model "${STAGE}/models/${model_id}.h5" \
      --fasta "${INPUT}/fixture/fixture.fa" \
      --windows "${STAGE}/fixture/mini_ccre_windows.tsv" \
      --fold "${INPUT}/prepared/actual_fold.json" \
      --output "${STAGE}/predictions/${model_id}" \
      --lineage hepatocyte \
      --batch-size 4 \
      --seed 20260824 \
      --expected-windows-per-role 4 \
      > "${STAGE}/validation/${model_id}.prediction.stdout.txt" \
      2> "${STAGE}/validation/${model_id}.prediction.stderr.txt"
done
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
  "${IMAGE}" python - "${CONFIG}" "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

import h5py
import numpy as np
import tensorflow as tf

from chrombpnet.training.utils.losses import multinomial_nll

config = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
root = Path(sys.argv[2]).resolve(strict=True)
budget = config["parameter_budget"]
parameter_counts = {}
for model_id in ("sequence_cnn_control", "sequence_transformer_control"):
    model = tf.keras.models.load_model(
        root / "models" / f"{model_id}.h5",
        custom_objects={"multinomial_nll": multinomial_nll, "tf": tf},
        compile=False,
    )
    if model.name != model_id or model.input_shape != (None, 2114, 4):
        raise SystemExit("trained model identity or input geometry differs")
    if [tuple(value) for value in model.output_shape] != [(None, 1000), (None, 1)]:
        raise SystemExit("trained model output geometry differs")
    parameter_counts[model_id] = model.count_params()
    if not budget["minimum_trainable_parameters"] <= model.count_params() <= budget["maximum_trainable_parameters"]:
        raise SystemExit("trained model parameter budget differs")
    receipt = json.loads(
        (root / "models" / f"{model_id}.training_receipt.json").read_text()
    )
    if receipt["test_outcomes_used"] or receipt["held_donor_atac_exposed"]:
        raise SystemExit("training receipt exposes held outcomes")
    summary = json.loads(
        (root / "predictions" / model_id / "summary.json").read_text()
    )
    contract = json.loads(
        (root / "predictions" / model_id / "prediction_contract.json").read_text()
    )
    if (
        summary["role_counts"] != {"valid": 4, "test": 4}
        or summary["observed_atac_input_exposed"]
        or summary["benchmark_metrics_calculated"]
        or not contract["ready_for_independent_evaluator_expansion"]
    ):
        raise SystemExit("fixed cCRE prediction contract differs")
    with h5py.File(
        root / "predictions" / model_id / "profile_probabilities.h5", "r"
    ) as handle:
        profiles = handle["profile_probability"][:]
    if profiles.shape != (8, 1000) or not np.allclose(
        profiles.sum(axis=1), 1.0, rtol=0, atol=1e-6
    ):
        raise SystemExit("fixed cCRE profile artifact differs")
    with (
        root / "predictions" / model_id / "standardized_predictions.tsv"
    ).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 8 or any(float(row["predicted"]) < 0 for row in rows):
        raise SystemExit("standardized prediction basis differs")
if max(parameter_counts.values()) / min(parameter_counts.values()) > budget["maximum_pairwise_ratio"]:
    raise SystemExit("control parameter-count parity differs")
validation = {
    "schema_version": "masld-bench-sequence-control-training-inference-probe-v1",
    "status": "pass",
    "models": ["sequence_cnn_control", "sequence_transformer_control"],
    "parameter_counts": parameter_counts,
    "seed": 20260824,
    "epochs": 1,
    "input_length": 2114,
    "output_length": 1000,
    "reverse_complement_tta": True,
    "held_donor_atac_exposed": False,
    "test_outcomes_used": False,
    "benchmark_metrics_calculated": False,
    "standardized_prediction_basis_passed": True,
}
(root / "validation" / "summary.json").write_text(
    json.dumps(validation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(validation, sort_keys=True))
PY
sha256sum \
  "${ROOT}/config/sequence_control_architectures.json" \
  "${ROOT}/scripts/sequence_control_architectures.py" \
  "${ROOT}/scripts/sequence_control_run_training.py" \
  "${ROOT}/scripts/sequence_control_predict_ccre.py" \
  "${ROOT}/slurm/sequence_control_probe_training_inference.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${INPUT}" \
  "${STAGE}/validation/input-verification.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

Path(sys.argv[2]).write_text(
    json.dumps(verify_frozen_tree(Path(sys.argv[1])), sort_keys=True, indent=2) + "\n",
    encoding="utf-8",
)
PY
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "sequence_control_training_inference_probe",
            "models": ["sequence_cnn_control", "sequence_transformer_control"],
            "input_artifacts_sha256": "5713097b91248632eb8486cd7e6c59d4c1f53832f69f5febca1631d113b480d2",
            "seed": 20260824,
            "epochs": 1,
            "held_donor_atac_exposed": False,
            "test_outcomes_used": False,
            "status": "passed",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Sequence-control training/inference probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
