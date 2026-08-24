#!/usr/bin/env bash
#SBATCH --job-name=masld-bpnet
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
case "${SEED}" in 20260824|20260825|20260826) ;;
  *) printf 'Invalid screening SEED\n' >&2; exit 2 ;;
esac
export PYTHONHASHSEED="${SEED}"
SPLIT_ID="donor${OUTER_FOLD}_genomic${OUTER_FOLD}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PSEUDOBULK="${ROOT}/executions/gse296875-training-fold-pseudobulk-21063993"
PSEUDOBULK_ARTIFACTS_SHA256=3190aa76cd5872f50d7b0c71bcb92520f1f2428110ce2fb4c26856994bc17e5b
INPUT="${ROOT}/executions/chrombpnet-${LINEAGE}-${SPLIT_ID}-inputs-${INPUT_JOB_ID}"
if [[ "${LINEAGE}:${OUTER_FOLD}:${INPUT_JOB_ID}" == hepatocyte:0:21064244 ]]; then
  INPUT="${ROOT}/executions/chrombpnet-hepatocyte-donor0-genomic0-inputs-21064244"
fi
if [[ "${INPUT_ARTIFACTS_SHA256}" == AUTO ]]; then
  if [[ ! -f "${INPUT}/ARTIFACTS.json" ]]; then
    printf 'Dependency-complete input receipt is missing: %s\n' "${INPUT}/ARTIFACTS.json" >&2
    exit 1
  fi
  read -r INPUT_ARTIFACTS_SHA256 _ < <(sha256sum "${INPUT}/ARTIFACTS.json")
fi
SPLIT="${ROOT}/executions/sequence-split-contract-v2-21063827"
REFERENCE="${ROOT}/executions/sequence-reference-contract-21063670"
CHROM_SIZES="${REFERENCE}/GRCh38.p14.sequence_model.chrom.sizes"
FASTA_ROOT="${ROOT}/executions/sequence-fasta-21063808"
FASTA="${FASTA_ROOT}/GRCh38.p14.sequence_model.fa"
BASE="${ROOT}/executions/environments/chrombpnet-native-tf2.8.2-21063230"
IMAGE="${BASE}/tensorflow-tensorflow-2.8.2-gpu.sif"
CORE="${ROOT}/executions/environments/chrombpnet-native-core-overlay-21063507"
ACQUISITION="${ROOT}/executions/upstream-acquisitions/bpnet-v2.0.0-21064186"
SOURCE="${ACQUISITION}/source/bpnet-f4593eedca51741f25b8c5cc0c0d647faf8c8a3c"
RUNTIME_PROBE="${ROOT}/executions/bpnet-native-runtime-probe-21064198"
TRAINING_PROBE="${ROOT}/executions/bpnet-training-inference-probe-21064291"
RESUME_PROBE="${ROOT}/executions/bpnet-training-resume-probe-21064293"
RESUME_PROBE_ARTIFACTS_SHA256=6f020b8df24764fd86c8c601697c0c5d40443cbe0c915e0f295a2652d83c3868
BIGWIG="${PSEUDOBULK}/donor_test_fold_${OUTER_FOLD}/${LINEAGE}.tn5.bw"
REGIONS="${INPUT}/regions/chrombpnet.regions.bed"
NONPEAKS="${INPUT}/nonpeaks/chrombpnet_negatives.bed"
FOLD="${INPUT}/regions/chrombpnet.fold.json"
WINDOWS="${SPLIT}/ccre_evaluation_windows.tsv"
TARGET="${ROOT}/executions/bpnet-${LINEAGE}-${SPLIT_ID}-seed${SEED}-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
PYTHONPATH_VALUE="${CORE}/site-packages:${SOURCE}:${ROOT}/scripts"
BPNET_CONTRACT_SHA256=0d806c472e70a31693b139fd0f065daf4416457f8b4d651a2a99e5aea0ff21b6
BPNET_TRAINING_SHA256=ca15d5ce83fc5d4f9e33c927e2ae9c1e63d6082db57aadcd91c1ea177b4765c8
BPNET_PREDICTION_SHA256=0e699cbee19598538af658ad3d26c3095035baa45cca264f2ec3af2b6deb8833
BPNET_ACQUISITION_VERIFY_SHA256=c93a6d3ab365989cbb4ae417f18216fc98c17c8d35959b38ef31cc892428f5bb

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed BPNet campaign: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite BPNet production campaign: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

printf '%s  %s\n' "${PSEUDOBULK_ARTIFACTS_SHA256}" \
  "${PSEUDOBULK}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${INPUT_ARTIFACTS_SHA256}" \
  "${INPUT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f \
  "${SPLIT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d \
  "${REFERENCE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 1901777f65ef618569eabb62967c5ace55a63b4bef8795af7a8838fc7042858b \
  "${CHROM_SIZES}" | sha256sum --check --strict
printf '%s  %s\n' 7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 \
  "${FASTA_ROOT}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' ed20dceeb315cffa2bc3de050d846eceefc552451fef741f0dd95b0149d10fce \
  "${FASTA}" | sha256sum --check --strict
printf '%s  %s\n' b52a78d42f568a260008b8e0bb541b2dc99aa816b844cf0b66e7551627f0a229 \
  "${BASE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 96495917c2d7e1b25d6f29b43f006c02d601314e50f250341e1f569fc51441aa \
  "${IMAGE}" | sha256sum --check --strict
printf '%s  %s\n' 61205e2ebbc012ddbeecbe957aba55d269e06c15d965c1e812a996c133a9f6f1 \
  "${CORE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' 8062574bfaa9067b0aca82ffa3eeae87188140ee5f5842399171ad4535cf88f6 \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' cd25afab907693c048eb71b0764db81984c5c664d0b267999682dc10cb3b05e4 \
  "${RUNTIME_PROBE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' ed43000954718d58870b55fcfed6ac9ae1fb0df1a853d777ac83c99cdb0ab051 \
  "${TRAINING_PROBE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${RESUME_PROBE_ARTIFACTS_SHA256}" \
  "${RESUME_PROBE}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${BPNET_CONTRACT_SHA256}" \
  "${ROOT}/scripts/bpnet_contract.py" | sha256sum --check --strict
printf '%s  %s\n' "${BPNET_TRAINING_SHA256}" \
  "${ROOT}/scripts/bpnet_run_training.py" | sha256sum --check --strict
printf '%s  %s\n' "${BPNET_PREDICTION_SHA256}" \
  "${ROOT}/scripts/bpnet_predict_ccre_windows.py" | sha256sum --check --strict
printf '%s  %s\n' "${BPNET_ACQUISITION_VERIFY_SHA256}" \
  "${ROOT}/scripts/bpnet_verify_acquisition.py" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${PSEUDOBULK}" "${INPUT}" "${SPLIT}" "${REFERENCE}" "${FASTA_ROOT}" \
  "${RUNTIME_PROBE}" "${TRAINING_PROBE}" "${RESUME_PROBE}" \
  "${STAGE}/validation/input-verification.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree

payload = {
    Path(value).name: verify_frozen_tree(Path(value)) for value in sys.argv[1:-1]
}
Path(sys.argv[-1]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY
python "${ROOT}/scripts/bpnet_verify_acquisition.py" \
  --acquisition "${ACQUISITION}" \
  --output "${STAGE}/validation/acquisition-verification.json"
python - "${INPUT}/validation/summary.json" \
  "${PSEUDOBULK}/training_pseudobulk_manifest.tsv" \
  "${OUTER_FOLD}" "${LINEAGE}" "${SPLIT_ID}" <<'PY'
from pathlib import Path
import csv
import json
import sys

summary = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
outer_fold = int(sys.argv[3])
lineage = sys.argv[4]
split_id = sys.argv[5]
if (
    summary["dataset_id"] != "gse296875"
    or summary["lineage_id"] != lineage
    or summary["split_id"] != split_id
    or summary["test_outcomes_available_to_training"]
    or summary["held_donor_validation_used_for_peak_discovery"]
    or summary["held_donor_test_used_for_peak_discovery"]
):
    raise SystemExit("shared training-input leakage contract differs")
with Path(sys.argv[2]).open(encoding="utf-8", newline="") as handle:
    rows = [
        row for row in csv.DictReader(handle, delimiter="\t")
        if int(row["donor_test_fold"]) == outer_fold and row["lineage_id"] == lineage
    ]
if (
    len(rows) != 1
    or int(rows[0]["donor_valid_fold"]) != (outer_fold + 1) % 5
    or len({int(value) for value in rows[0]["donor_train_folds"].split(",")}) != 3
    or int(rows[0]["donors"]) < 1
):
    raise SystemExit("BPNet training pseudobulk manifest differs")
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
  --env OMP_NUM_THREADS=4 \
  --env OPENBLAS_NUM_THREADS=4 \
  --env MKL_NUM_THREADS=4 \
  "${IMAGE}" python "${ROOT}/scripts/bpnet_run_training.py" \
    --genome "${FASTA}" \
    --chrom-sizes "${CHROM_SIZES}" \
    --bigwig "${BIGWIG}" \
    --regions "${REGIONS}" \
    --nonpeaks "${NONPEAKS}" \
    --fold "${FOLD}" \
    --output "${STAGE}/model" \
    --epochs 50 \
    --batch-size 64 \
    --threads 4 \
    --learning-rate 0.001 \
    --min-learning-rate 0.00001 \
    --max-jitter 128 \
    --negative-ratio 0.1 \
    --early-stop 5 \
    --early-stop-min-delta 0.001 \
    --reduce-lr-patience 2 \
    --reduce-lr-factor 0.5 \
    --seed "${SEED}" \
    > "${STAGE}/validation/training.stdout.txt" \
    2> "${STAGE}/validation/training.stderr.txt"

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  --env PYTHONDONTWRITEBYTECODE=1 \
  --env "PYTHONHASHSEED=${SEED}" \
  --env TF_FORCE_GPU_ALLOW_GROWTH=true \
  "${IMAGE}" python "${ROOT}/scripts/bpnet_predict_ccre_windows.py" \
    --model "${STAGE}/model/bpnet.inference.h5" \
    --fasta "${FASTA}" \
    --windows "${WINDOWS}" \
    --fold "${FOLD}" \
    --output "${STAGE}/predictions" \
    --batch-size 256 \
    --expected-windows 32000 \
    --seed "${SEED}" \
    > "${STAGE}/validation/inference.stdout.txt" \
    2> "${STAGE}/validation/inference.stderr.txt"

singularity exec --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${CORE}/site-packages:${ROOT}/scripts" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python - "${STAGE}" "${LINEAGE}" "${SPLIT_ID}" "${SEED}" <<'PY'
from pathlib import Path
import csv
import json
import sys

import h5py
import numpy as np
import tensorflow as tf

root = Path(sys.argv[1]).resolve(strict=True)
lineage = sys.argv[2]
split_id = sys.argv[3]
seed = int(sys.argv[4])
training = json.loads((root / "model" / "training_receipt.json").read_text())
prediction = json.loads((root / "predictions" / "summary.json").read_text())
if (
    not training["campaign_complete"]
    or training["requested_epochs"] != 50
    or training["architecture"]["parameter_count"] != 109124
    or training["architecture"]["bias_or_control_input"]
    or training["test_contigs_exposed"]
    or training["observed_held_donor_atac_exposed"]
    or training["evaluator_outcomes_exposed"]
    or training["max_jitter"] != 128
    or training["boundary_policy"] != {
        "training_required_flank_bp": 1185,
        "validation_required_flank_bp": 1057,
        "input_half_width_bp": 1057,
        "max_jitter_bp": 128,
        "role_specific": True,
        "applied_to_peaks_and_nonpeaks": True,
    }
):
    raise SystemExit("BPNet training receipt differs")
if prediction["role_counts"] != {"valid": 16000, "test": 16000}:
    raise SystemExit("BPNet prediction role census differs")
model = tf.keras.models.load_model(root / "model" / "bpnet.inference.h5", compile=False)
if model.input_shape != (None, 2114, 4) or len(model.inputs) != 1:
    raise SystemExit("BPNet exported input contract differs")
with h5py.File(root / "predictions" / "profile_probabilities.h5") as handle:
    profiles = handle["profile_probability"]
    if profiles.shape != (32000, 1000):
        raise SystemExit("BPNet profile prediction shape differs")
    for start in range(0, profiles.shape[0], 1024):
        values = profiles[start : start + 1024]
        if not np.isfinite(values).all():
            raise SystemExit("BPNet profile prediction is non-finite")
        np.testing.assert_allclose(values.sum(axis=1), 1.0, rtol=0, atol=1e-6)
with (root / "predictions" / "regional_counts.tsv").open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if len(rows) != 32000 or any(
    float(row["strand_averaged_mass"]) < 0 for row in rows
):
    raise SystemExit("BPNet regional count output differs")
(root / "validation" / "summary.json").write_text(
    json.dumps({
        "schema_version": "masld-bench-bpnet-production-campaign-v1",
        "status": "pass",
        "model_id": "bpnet",
        "dataset_id": "gse296875",
        "lineage_id": lineage,
        "split_id": split_id,
        "seed": seed,
        "completed_epochs": training["completed_epochs"],
        "best_epoch": training["best_epoch"],
        "parameter_count": training["architecture"]["parameter_count"],
        "validation_prediction_windows": 16000,
        "test_prediction_windows": 16000,
        "bias_or_control_input": False,
        "held_donor_outcomes_used_for_training": False,
        "evaluator_outcomes_exposed": False,
        "benchmark_metrics_calculated": False,
    }, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY

singularity exec --nv --cleanenv \
  --bind "${ROOT}:${ROOT}:rw" \
  --env "PYTHONPATH=${PYTHONPATH_VALUE}" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/bpnet_contract.py" \
  "${ROOT}/scripts/bpnet_run_training.py" \
  "${ROOT}/scripts/bpnet_predict_ccre_windows.py" \
  "${ROOT}/slurm/bpnet_train_predict_hepatocyte_donor0_genomic0.sbatch" \
  > "${STAGE}/validation/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" \
  "${PSEUDOBULK_ARTIFACTS_SHA256}" "${INPUT_ARTIFACTS_SHA256}" \
  "${LINEAGE}" "${SPLIT_ID}" "${SEED}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(freeze_tree(root, {
    "artifact_class": "bpnet_training_prediction_campaign",
    "model_id": "bpnet",
    "dataset_id": "gse296875",
    "lineage_id": sys.argv[4],
    "split_id": sys.argv[5],
    "seed": int(sys.argv[6]),
    "source_revision": "f4593eedca51741f25b8c5cc0c0d647faf8c8a3c",
    "pseudobulk_artifacts_sha256": sys.argv[2],
    "shared_training_input_artifacts_sha256": sys.argv[3],
    "split_artifacts_sha256": "00c073e4667c16a4c54dc013c57d2812dbc88b8ecacd5cb8a6088fca9072fe0f",
    "reference_artifacts_sha256": "7e761c93a6e6aafcc54decb611a8462d2ea5a7c9466e41f153173acee6d8b29d",
    "fasta_artifacts_sha256": "7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39",
    "runtime_probe_artifacts_sha256": "cd25afab907693c048eb71b0764db81984c5c664d0b267999682dc10cb3b05e4",
    "bias_or_control_input": False,
    "evaluator_outcomes_exposed": False,
    "benchmark_metrics_calculated": False,
    "status": "passed",
}))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'BPNet training/prediction campaign: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
