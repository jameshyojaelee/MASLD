#!/usr/bin/env bash
#SBATCH --job-name=model-work-335
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
: "${MASLD_GPU_DISPATCHER_AUTHORITY:?dispatcher authority is required}"
: "${QUERY_FIXTURE_ROOT:?query fixture root is required}"
: "${QUERY_FIXTURE_ARTIFACTS_SHA256:?query fixture hash is required}"
if [[ ${MASLD_GPU_DISPATCHER_AUTHORITY} != sequence_gpu_dispatcher_v2 ]]; then
  printf 'Invalid GPU dispatcher authority.\n' >&2
  exit 1
fi

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
CHECKPOINT=${ROOT}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens
SCRIPT=${ROOT}/scripts/transcriptformer_extract_unlabeled_embeddings.py
TARGET=${ROOT}/executions/model-work-335-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  status=$?
  if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
  exit "${status}"
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite TF-Sapiens query embeddings.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in \
  "${RUNTIME_SHA256} ${RUNTIME}/ARTIFACTS.json" \
  "eff027e7393f92a32cadcbb2e13a917e24b62b71341571bf40f38ff6ec6d3021 ${CHECKPOINT}/model_weights.pt" \
  "2413d512e8fdd550389bc574dc465357e943e23b4411cab6d23423af08095c2f ${CHECKPOINT}/config.json" \
  "676640be87ba4ba7ea7af8d9dc470a8ae8be02768989478027dec5fe61f465e2 ${CHECKPOINT}/vocabs/homo_sapiens_gene.h5" \
  "0c405e4ead45a4b8350d8e874f834273efdd3f7ad4669b52d3e3727fb4fe70af ${CHECKPOINT}/vocabs/assay_vocab.json" \
  "${QUERY_FIXTURE_ARTIFACTS_SHA256} ${QUERY_FIXTURE_ROOT}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${RUNTIME}/env/bin/python" - <<'PY' "${QUERY_FIXTURE_ROOT}"
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
meta = verify_frozen_tree(Path(sys.argv[1]))["metadata"]
if meta.get("artifact_class") != "gse296875_transcriptformer_outcome_blind_query_fixture":
    raise SystemExit("query fixture authority differs")
if any(meta.get(key) is not False for key in ("labels_read", "donor_ids_read", "atac_read", "target_adaptation_performed", "sealed_outcomes_read")):
    raise SystemExit("query fixture firewall differs")
PY
"${RUNTIME}/env/bin/python" "${SCRIPT}" \
  --checkpoint-root "${CHECKPOINT}" \
  --fixture-root "${QUERY_FIXTURE_ROOT}" \
  --output "${STAGE}/embeddings" \
  --expected-rows 7500 \
  > "${STAGE}/validation/extraction.stdout.json" \
  2> "${STAGE}/validation/extraction.stderr.txt"
"${RUNTIME}/env/bin/python" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
nvidia-smi -q > "${STAGE}/validation/nvidia_smi.txt"
sha256sum "${SCRIPT}" slurm/run_gse296875_transcriptformer_query_embeddings.sbatch \
  "${QUERY_FIXTURE_ROOT}/ARTIFACTS.json" "${RUNTIME}/ARTIFACTS.json" \
  "${CHECKPOINT}/model_weights.pt" "${CHECKPOINT}/config.json" \
  "${CHECKPOINT}/vocabs/homo_sapiens_gene.h5" \
  "${CHECKPOINT}/vocabs/assay_vocab.json" > "${STAGE}/source.sha256"
"${RUNTIME}/env/bin/python" - <<'PY' "${STAGE}"
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1])
receipt = json.loads((root / "embeddings/receipt.json").read_text())
if receipt.get("rows") != 7500 or receipt.get("evaluation_label_columns_read") != [] or receipt.get("sealed_outcomes_read") is not False:
    raise SystemExit("query embedding receipt differs")
freeze_tree(root, {
    "artifact_class": "gse296875_transcriptformer_outcome_blind_embeddings",
    "view_id": "gse296875_rna_cell_state_external_development_7500_v1",
    "model_id": "transcriptformer_tf_sapiens",
    "rows": 7500,
    "labels_read": False,
    "donor_ids_read": False,
    "atac_read": False,
    "target_adaptation_performed": False,
    "metrics_calculated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
