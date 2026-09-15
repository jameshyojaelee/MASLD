#!/usr/bin/env bash
#SBATCH --job-name=model-data-335
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
FEATURES=${ROOT}/executions/model-data-332-21100842/view/features
FEATURES_SHA256=1cacb52032beb27969c69ae2d00cfdf3d18b333f1458bb4fb830f5225d79dc6c
VOCAB=${ROOT}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens/vocabs/homo_sapiens_gene.h5
VOCAB_SHA256=676640be87ba4ba7ea7af8d9dc470a8ae8be02768989478027dec5fe61f465e2
TARGET=${ROOT}/executions/model-data-335-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
PYTHON=${RUNTIME}/env/bin/python

preserve_failed_stage() {
  status=$?
  if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
  exit "${status}"
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite query-fixture campaign.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in \
  "${RUNTIME_SHA256} ${RUNTIME}/ARTIFACTS.json" \
  "${FEATURES_SHA256} ${FEATURES}/ARTIFACTS.json" \
  "${VOCAB_SHA256} ${VOCAB}"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
"${PYTHON}" scripts/build_gse296875_transcriptformer_query_fixture.py \
  --features "${FEATURES}" --expected-features-sha256 "${FEATURES_SHA256}" \
  --vocab "${VOCAB}" --output "${STAGE}/fixture" \
  > "${STAGE}/validation/build.stdout.json" \
  2> "${STAGE}/validation/build.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum scripts/build_gse296875_transcriptformer_query_fixture.py \
  slurm/build_gse296875_transcriptformer_query_fixture_cpu.sbatch \
  "${FEATURES}/ARTIFACTS.json" "${RUNTIME}/ARTIFACTS.json" "${VOCAB}" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - <<'PY' "${STAGE}"
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1])
fixture = root / "fixture"
meta = verify_frozen_tree(fixture)["metadata"]
receipt = json.loads((fixture / "fixture_receipt.json").read_text())
if (
    meta.get("artifact_class") != "gse296875_transcriptformer_outcome_blind_query_fixture"
    or meta.get("rows") != 7500
    or receipt.get("policies", {}).keys() != {"native", "common"}
    or any(meta.get(key) is not False for key in ("labels_read", "donor_ids_read", "atac_read", "target_adaptation_performed", "sealed_outcomes_read"))
):
    raise SystemExit("query fixture campaign differs")
freeze_tree(root, {
    "artifact_class": "gse296875_transcriptformer_query_fixture_campaign",
    "view_id": "gse296875_rna_cell_state_external_development_7500_v1",
    "rows": 7500,
    "labels_read": False,
    "donor_ids_read": False,
    "atac_read": False,
    "target_adaptation_performed": False,
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/fixture/ARTIFACTS.json"
