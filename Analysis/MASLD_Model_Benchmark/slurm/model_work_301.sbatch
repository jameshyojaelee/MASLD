#!/usr/bin/env bash
#SBATCH --job-name=model-work-301
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=20:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
if [[ "${MASLD_GPU_DISPATCHER_AUTHORITY}" != sequence_gpu_dispatcher_v2 ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 1
fi
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/caduceus-native-torch2.2.0-21065466"
CHECKPOINTS="${ROOT}/executions/dna-lm-open-safetensors-21064556"
CODE="${ROOT}/executions/dna-lm-custom-code-21064626"
COMMON_PROBE="${ROOT}/executions/caduceus-native-probe-21065951"
LONG_PROBE="${ROOT}/executions/caduceus-131k-probe-21066018"
FIXTURE="${ROOT}/executions/gse281364-caduceus-131k-fixture-21071004"
FASTA="${ROOT}/executions/sequence-fasta-21063808"
PREFLIGHT="${ROOT}/executions/caduceus-131k-gse281364-preflight-21072369"
CONFIG="${ROOT}/config/gse281364_caduceus_131k_contract.json"
CACHE="${ROOT}/executions/resumable-model-work-301"
CACHE_ID=caduceus-131k-gse281364-v1
TARGET="${ROOT}/executions/model-work-301-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
RUNTIME_TMP=$(mktemp -d "${TMPDIR:-/tmp}/model-work-301-${SLURM_JOB_ID}.XXXXXX")
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
  if [[ -d "${RUNTIME_TMP}" ]]; then rm -r "${RUNTIME_TMP}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite model-work-301 artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for binding in \
  "3d59792d8a13ca10a2df3d6386bc266eb0e108def6279188793d1cfee4fc1167 ${RUNTIME}" \
  "6d21df999835624905e08573038dc72e893a5b25a8ca8f4ec1acea8f61d91c77 ${CHECKPOINTS}" \
  "70906faae158909a8972f674f680943d8a51467032ae3266a6a12f2152beb36c ${CODE}" \
  "6bc5288848ca6cd696cc90c36ceea4e789e97b7f11ceb70670526bd0ee696ceb ${COMMON_PROBE}" \
  "a9c531b8a044c8e7990e24cbafcb7b035d51d5861efccb7c2c326c66f41b5456 ${LONG_PROBE}" \
  "82b47ba2059193cae23e35db5401013659879024052c47558e602c03b4ee2c70 ${FIXTURE}" \
  "7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 ${FASTA}" \
  "08990558a6abfae5cf57f21b6a97f792c21ecef443717061e581908be024d562 ${PREFLIGHT}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" | sha256sum --check --strict
done
for binding in \
  "0c3e50fb5dc5a41c3686caf638416106b523a1e059a28c50fd5a75c2e37a0dbf ${CONFIG}" \
  "e3d2cb3b7ccc209feb77d3cbd75e6d68e270be96d5f3cfcfef5de1dab1b618de ${ROOT}/scripts/caduceus_predict_gse281364_131k.py" \
  "6405c4dbc188e8d904b8e5d5afba68801e824773ac06becd462dd88b4327821f ${ROOT}/scripts/project_gse281364_caduceus_131k_embeddings.py" \
  "3aa9c27842c351af29d68f7e8642815e442eca26a1b470ccadffd53dd2d5876d ${ROOT}/tests/unit/test_caduceus_predict_gse281364_131k.py"; do
  read -r expected path <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
install -d -m 0750 "${STAGE}/raw" "${STAGE}/validation" \
  "${STAGE}/offline-guard" "${STAGE}/hf-home" "${CACHE}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${CHECKPOINTS}" "${CODE}" \
  "${COMMON_PROBE}" "${LONG_PROBE}" "${FIXTURE}" "${FASTA}" "${PREFLIGHT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
module purge
cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket
_socket = socket.socket
class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for model-work-301")
        return super().__new__(cls, family, *args, **kwargs)
socket.socket = GuardedSocket
def blocked_connection(*args, **kwargs):
    raise RuntimeError("network disabled for model-work-301")
socket.create_connection = blocked_connection
PY
tar -xzf "${RUNTIME}/payload/caduceus-native-env.tar.gz" -C "${RUNTIME_TMP}"
export HF_HOME="${STAGE}/hf-home"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export PYTHONPATH="${STAGE}/offline-guard:${ROOT}"
cd "${ROOT}"
"${RUNTIME_TMP}/bin/python" -m unittest -v \
  tests.unit.test_caduceus_predict_gse281364_131k \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${RUNTIME_TMP}/bin/python" scripts/caduceus_predict_gse281364_131k.py \
  --config "${CONFIG}" \
  --code-root "${CODE}/code/caduceus" \
  --checkpoint "${CHECKPOINTS}/checkpoints/caduceus.model.safetensors" \
  --fixture-artifacts "${FIXTURE}/ARTIFACTS.json" \
  --manifest "${FIXTURE}/fixture/sequence_manifest.tsv" \
  --fasta-artifacts "${FASTA}/ARTIFACTS.json" \
  --fasta "${FASTA}/GRCh38.p14.sequence_model.fa" \
  --fai "${FASTA}/GRCh38.p14.sequence_model.fa.fai" \
  --checkpoint-root "${CACHE}" \
  --cache-id "${CACHE_ID}" \
  --chunk-size 8 \
  --output "${STAGE}/raw/caduceus" \
  > "${STAGE}/validation/prediction.stdout.jsonl" \
  2> "${STAGE}/validation/prediction.stderr.txt"
"${RUNTIME_TMP}/bin/python" scripts/project_gse281364_caduceus_131k_embeddings.py \
  --raw-root "${STAGE}/raw/caduceus" \
  --output "${STAGE}/projected" \
  > "${STAGE}/validation/projection.stdout.json" \
  2> "${STAGE}/validation/projection.stderr.txt"
"${RUNTIME_TMP}/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
nvidia-smi -q > "${STAGE}/validation/nvidia-smi-q.txt"
module list 2> "${STAGE}/validation/module-list.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/caduceus_predict_gse281364_131k.py" \
  "${ROOT}/scripts/project_gse281364_caduceus_131k_embeddings.py" \
  "${ROOT}/tests/unit/test_caduceus_predict_gse281364_131k.py" \
  "${ROOT}/slurm/model_work_301.sbatch" \
  "${PREFLIGHT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
rm -r "${RUNTIME_TMP}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
raw = json.loads((root / "raw/caduceus/receipt.json").read_text())
projection = json.loads((root / "projected/receipt.json").read_text())
tests = (root / "validation/tests.stderr.txt").read_text()
if (
    raw["status"] != "pass_outcome_blind_embedding_extraction"
    or raw["elements"] != 1033
    or raw["outer_locus_sequence_groups"] != 1033
    or raw["outer_folds"] != 5
    or raw["input_length_bp"] != 131072
    or raw["hidden_width"] != 256
    or raw["max_external_rc_feature_abs_difference"] > 1e-5
    or raw["deterministic_repeat_max_abs_difference"] > 1e-5
    or raw["license"] != "Apache-2.0"
    or raw["exposure_state"] != "target_label_unexposed"
    or raw["restricted_comparator"]
    or not raw["open_champion_eligible_after_task_and_external_evaluation_gates"]
    or raw["reporter_counts_read"]
    or raw["reporter_outcomes_read"]
    or raw["sealed_labels_read"]
    or raw["downstream_head_fit"]
    or raw["runtime_network_allowed"]
    or projection["status"] != "pass_outcome_blind_fold_safe_projection"
    or projection["elements"] != 1033
    or projection["projection_width"] != 256
    or projection["head_input_width"] != 1024
    or projection["projection_fit_on_held_out_fold"]
    or projection["reporter_outcomes_read"]
    or projection["sealed_labels_read"]
    or projection["downstream_head_fit"]
    or "Ran 3 tests" not in tests
    or "OK" not in tests
):
    raise SystemExit("model-work-301 receipt differs")
freeze_tree(root, {
    "artifact_class": "outcome_blind_caduceus_131k_embeddings",
    "dataset_id": "gse281364",
    "model_id": "caduceus",
    "checkpoint_sha256": "a3e6976fe90460ff5d90d371b457d0263dc65e156ea5651b4a452e98228daef2",
    "license": "Apache-2.0",
    "exposure_state": "target_label_unexposed",
    "elements": 1033,
    "outer_locus_sequence_groups": 1033,
    "outer_folds": 5,
    "input_length_bp": 131072,
    "projection_width": 256,
    "head_input_width": 1024,
    "reporter_outcomes_read": False,
    "downstream_head_fit": False,
    "champion_eligible_only_after_task_and_external_evaluation_gates": True,
    "status": "passed",
})
PY
module purge
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'model-work-301 artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
