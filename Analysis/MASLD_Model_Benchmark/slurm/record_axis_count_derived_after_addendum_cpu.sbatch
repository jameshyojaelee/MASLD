#!/usr/bin/env bash
#SBATCH --job-name=model-check-1048
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 0b, job 3 of 3: two measurements the frozen check did not make.
#
# Additive overlay. It alters no frozen verdict and registers no check. Both
# measurements were commissioned AFTER model-check-1047 had already computed and
# held-back the outcome, so they are recorded as derived-after and may not be cited
# as pre-registered. The frozen result tree is pinned by digest below so the
# ordering is provable from the output files.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/record_axis_count_derived_after_addendum.py
ANALYSIS=${ROOT}/scripts/evaluate_gse267145_axis_count.py
RESULT=${ROOT}/executions/model-check-1047-21181822-axis-count/axis_count
ENDPOINTS=${ROOT}/executions/model-data-061-21079623/activation/participant_endpoints.tsv
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular

TARGET=${ROOT}/executions/model-check-1048-${SLURM_JOB_ID}-axis-count-derived-after
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the addendum record.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/record_axis_count_derived_after_addendum_cpu.sbatch

cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
ae044ef6285498824105acf85945583e6d44f4c60445cc39bcf811bea7269afc  executions/model-check-1047-21181822-axis-count/ARTIFACTS.json
78ce955b6ef43ff45ebf78ec7744aa464febf3da93ce8aff4b6daf0f80c4ea3e  executions/model-check-1047-21181822-axis-count/axis_count/ARTIFACTS.json
e7397c3038dfbbdd1619048acfdd6cabc0e018a7af161175925130a7a6ffed3e  executions/model-check-1046-21181811-axis-count-prespec/prespecification/ARTIFACTS.json
9a13c0fd65db0483566b7495def11ade3339861de091ea29fcbdc4c34f47b5b9  executions/model-data-061-21079623/activation/participant_endpoints.tsv
812d5fce6e7cc88f6939ae493db7b7e7470057c21cf3b219abc39f587ff64ddc  executions/model-data-064-21079902/fixture/molecular/rna_values.npy
baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd  executions/model-data-064-21079902/fixture/molecular/rna_feature_axis.tsv
7e5be1a16df505e5daad451398e8b1a4b3cbbc48a32cd1381113e1ae4b7347bc  executions/model-data-064-21079902/fixture/molecular/participant_axis.tsv
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" - "${RESULT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
print("frozen Stage 0b result verified")
PY

"${PYTHON}" -m unittest -v tests.unit.test_record_axis_count_derived_after_addendum \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" "${SCRIPT}" \
  --analysis-script "${ANALYSIS}" \
  --result "${RESULT}" \
  --endpoints "${ENDPOINTS}" \
  --molecular "${MOLECULAR}" \
  --output "${STAGE}/addendum" \
  --count-null-draws 20000 \
  --null-block 250 \
  --splits 200 \
  --seed 20260827 \
  > "${STAGE}/addendum_stdout.json" 2> "${STAGE}/addendum_stderr.txt"

sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" "${ANALYSIS}" \
  slurm/record_axis_count_derived_after_addendum_cpu.sbatch \
  tests/unit/test_record_axis_count_derived_after_addendum.py \
  "${RESULT}/axis_count.json" "${ENDPOINTS}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
payload = json.loads(
    (root / "addendum" / "derived_after_addendum.json").read_text(encoding="utf-8")
)
for key in ("record_is_additive_overlay",):
    if payload[key] is not True:
        raise SystemExit(f"{key} is not asserted")
for key in ("frozen_verdict_altered", "in_force_for_the_frozen_decision", "is_a_gate"):
    if payload[key] is not False:
        raise SystemExit(f"{key} must be false for an overlay")
ordering = payload["honesty_about_ordering"]
if ordering["frozen_before_the_analysis"] is not False:
    raise SystemExit("the overlay claims to predate the analysis")
if ordering["recorded_after_the_analysis"] is not True:
    raise SystemExit("the overlay does not disclose that it came after")
if ordering["may_not_be_cited_as_pre_registered"] is not True:
    raise SystemExit("the overlay does not forbid citing it as pre-registered")
if payload["saf_activity_versus_fibrosis_partial"]["is_a_gate"] is not False:
    raise SystemExit("the A-vs-F partial is marked a gate")
if payload["frozen_artifacts_this_record_covers"]["outcome_in_the_frozen_result"] != (
    "ONE_AXIS_ACTIVITY"
):
    raise SystemExit("the covered outcome is not the frozen one")

verify_frozen_tree(root / "addendum")
freeze_tree(root, {
    "artifact_class": "gse267145_axis_count_derived_after_campaign",
    "stage": "stage_0b_gate",
    "job": "3_of_3",
    "record_is_additive_overlay": True,
    "frozen_verdict_altered": False,
    "is_a_gate": False,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "outcome_unchanged": payload["frozen_artifacts_this_record_covers"][
        "outcome_in_the_frozen_result"],
    "reliability_prediction_held": payload["split_half_reliability"]["prediction_held"],
    "saf_activity_reliability": payload["split_half_reliability"][
        "measured_saf_activity_full_length_reliability"],
}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/addendum/ARTIFACTS.json"
