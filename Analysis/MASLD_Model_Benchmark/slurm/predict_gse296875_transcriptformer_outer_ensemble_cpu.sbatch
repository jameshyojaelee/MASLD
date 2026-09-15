#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-transfer-337
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
QUERY=${ROOT}/executions/model-work-335-21101384
SOURCE=${ROOT}/executions/model-training-500-21080428
SOURCE_SHA256=fc4601988491b086d5e89f0ff35f2f73d9c27cf74846a6ff7895efc130c825a5
HEADS=${ROOT}/executions/model-training-21082720
HEADS_SHA256=f0cefd02290c2b2135c3578784076257e2a63213a99e670bf848a47f18a1e27e
SPLIT=${ROOT}/executions/model-data-052-21077595/split
SPLIT_SHA256=10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680
TARGET=${ROOT}/executions/model-cpu-transfer-337-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
trap 'status=$?; if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi; exit ${status}' EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then exit 1; fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in "${RUNTIME_SHA256} ${RUNTIME}/ARTIFACTS.json" "${SOURCE_SHA256} ${SOURCE}/ARTIFACTS.json" "${HEADS_SHA256} ${HEADS}/ARTIFACTS.json" "${SPLIT_SHA256} ${SPLIT}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
QUERY_SHA256=$(sha256sum "${QUERY}/ARTIFACTS.json" | cut -d ' ' -f 1)
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${RUNTIME}/env/bin/python" - <<'PY' "${QUERY}" "${SOURCE}" "${HEADS}" "${SPLIT}"
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]: verify_frozen_tree(Path(value))
PY
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${RUNTIME}/env/bin/python" scripts/predict_gse296875_transcriptformer_outer_ensemble.py \
  --query-embeddings "${QUERY}" --expected-query-embeddings-sha256 "${QUERY_SHA256}" \
  --source-embeddings "${SOURCE}" --expected-source-embeddings-sha256 "${SOURCE_SHA256}" \
  --source-heads "${HEADS}" --expected-source-heads-sha256 "${HEADS_SHA256}" \
  --split "${SPLIT}" \
  --output "${STAGE}/branch" \
  > "${STAGE}/validation/predict.stdout.json" 2> "${STAGE}/validation/predict.stderr.txt"
"${RUNTIME}/env/bin/python" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum scripts/predict_gse296875_transcriptformer_outer_ensemble.py \
  slurm/predict_gse296875_transcriptformer_outer_ensemble_cpu.sbatch \
  "${QUERY}/ARTIFACTS.json" "${SOURCE}/ARTIFACTS.json" "${HEADS}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${RUNTIME}/env/bin/python" - <<'PY' "${STAGE}"
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]); meta = verify_frozen_tree(root / "branch")["metadata"]
if (
    meta.get("branch_id") != "transcriptformer_outer_ensemble"
    or meta.get("models") != ["transcriptformer_tf_sapiens__linear", "transcriptformer_tf_sapiens__two_layer_mlp"]
    or meta.get("source_fit_census_per_model") != 15
):
    raise SystemExit("TF-Sapiens transfer branch differs")
freeze_tree(root, {
    "artifact_class": "gse296875_transcriptformer_transfer_campaign",
    "view_id": "gse296875_rna_cell_state_external_development_7500_v1",
    "models": meta["models"],
    "rows": 7500,
    "query_labels_read": False,
    "query_donor_ids_read": False,
    "query_atac_read": False,
    "metrics_calculated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/branch/ARTIFACTS.json"
