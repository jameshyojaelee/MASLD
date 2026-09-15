#!/usr/bin/env bash
#SBATCH --job-name=model-check-1047
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 0b of the MASLD showcase model, job 2 of 2: the axis-count check.
#
# The prespecification frozen by model-check-1046 (job 21181811) is verified by
# digest BEFORE a single expression value is read. The pinned digests below were
# read from that completed tree, so the ordering is provable from the output files
# and does not rest on a claim in a report. Stage 0's own prespecification and
# result trees are pinned too, because the frozen document cites them and a
# citation that is not checked is prose.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/evaluate_gse267145_axis_count.py

PRESPEC_ROOT=${ROOT}/executions/model-check-1046-21181811-axis-count-prespec
PRESPEC=${PRESPEC_ROOT}/prespecification
ENDPOINTS=${ROOT}/executions/model-data-061-21079623/activation/participant_endpoints.tsv
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular
FOLDS=${ROOT}/executions/model-data-064-21079902/fixture/folds/participant_outer_folds.tsv

TARGET=${ROOT}/executions/model-check-1047-${SLURM_JOB_ID}-axis-count
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the axis-count result.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/evaluate_gse267145_axis_count_cpu.sbatch

# ---------------------------------------------------------------------------
# The prespecification came first. Prove it before opening the matrix.
# ---------------------------------------------------------------------------
cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
217bdc325cc5205ce007857873e33c6c0dc13e626ae862d4de235acb90226fd5  executions/model-check-1046-21181811-axis-count-prespec/ARTIFACTS.json
e7397c3038dfbbdd1619048acfdd6cabc0e018a7af161175925130a7a6ffed3e  executions/model-check-1046-21181811-axis-count-prespec/prespecification/ARTIFACTS.json
330411c86ae47452efd293ab8c9620813ddfd55d1e774ea86e62006ba7be39ac  executions/model-check-1046-21181811-axis-count-prespec/prespecification/axis_count_prespecification.json
5395050a71a8e7a8bdc9ecd6b126c08e121a2f4ff6f42268bed0cbcdfe3efa76  executions/model-check-1043-21181265-aspect-separability-prespec/prespecification/ARTIFACTS.json
682238189f8bd97144eb7a7d39efefea8f17546d60d25ad655498af992f60c48  executions/model-check-1044-21181334-aspect-separability/separability/ARTIFACTS.json
9a13c0fd65db0483566b7495def11ade3339861de091ea29fcbdc4c34f47b5b9  executions/model-data-061-21079623/activation/participant_endpoints.tsv
812d5fce6e7cc88f6939ae493db7b7e7470057c21cf3b219abc39f587ff64ddc  executions/model-data-064-21079902/fixture/molecular/rna_values.npy
baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd  executions/model-data-064-21079902/fixture/molecular/rna_feature_axis.tsv
7e5be1a16df505e5daad451398e8b1a4b3cbbc48a32cd1381113e1ae4b7347bc  executions/model-data-064-21079902/fixture/molecular/participant_axis.tsv
9ddfaf9c58b58a5b239eebbbc19b214210543fb4f60f3f5a7c563a70fbebb0d1  executions/model-data-064-21079902/fixture/folds/participant_outer_folds.tsv
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

# The prespecification tree must also be internally consistent, not merely
# match a digest this wrapper carries.
"${PYTHON}" - "${PRESPEC}" "${PRESPEC_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for argument in sys.argv[1:]:
    verify_frozen_tree(Path(argument).resolve(strict=True))
print("prespecification tree verified")
PY

"${PYTHON}" -m unittest -v tests.unit.test_evaluate_gse267145_axis_count \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

find "${MOLECULAR}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_before.txt"

# The script owns its --output directory and refuses to overwrite one.
"${PYTHON}" "${SCRIPT}" \
  --prespecification "${PRESPEC}" \
  --endpoints "${ENDPOINTS}" \
  --molecular "${MOLECULAR}" \
  --folds "${FOLDS}" \
  --output "${STAGE}/axis_count" \
  --count-null-draws 20000 \
  --null-block 250 \
  --label-draws 20000 \
  --bootstrap 1000 \
  --seed 20260827 \
  > "${STAGE}/axis_count_stdout.txt" 2> "${STAGE}/axis_count_stderr.txt"

find "${MOLECULAR}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_after.txt"
cmp -s "${STAGE}/frozen_before.txt" "${STAGE}/frozen_after.txt" || { printf 'A frozen artifact changed.\n' >&2; exit 1; }
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" \
  slurm/evaluate_gse267145_axis_count_cpu.sbatch \
  tests/unit/test_evaluate_gse267145_axis_count.py \
  "${PRESPEC}/axis_count_prespecification.json" \
  "${ENDPOINTS}" "${FOLDS}" \
  "${MOLECULAR}/rna_values.npy" "${MOLECULAR}/rna_feature_axis.tsv" \
  "${MOLECULAR}/participant_axis.tsv" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
result = root / "axis_count"
payload = json.loads((result / "axis_count.json").read_text(encoding="utf-8"))

if payload["prespec_id"] != "gse267145_axis_count_v1":
    raise SystemExit("the result is not against the frozen prespecification")
declared = {
    "THREE_AXIS_SAF", "TWO_AXIS_ACTIVITY_FIBROSIS", "TWO_AXIS_STEATOSIS_ACTIVITY",
    "ONE_AXIS_ACTIVITY", "ONE_AXIS_FIBROSIS", "ONE_AXIS_UNDETERMINED",
    "INDETERMINATE",
}
if payload["outcome"] not in declared:
    raise SystemExit(f"undeclared outcome {payload['outcome']!r}")
if payload["a_one_axis_answer_is_a_result_not_a_failure"] is not True:
    raise SystemExit("a one-axis answer is not recorded as a result")
if payload["cohort"]["distinct_participants"] != 99:
    raise SystemExit("the participant count is not 99")
if payload["honest_provenance"][
    "this_question_was_derived_after_seeing_stage_0s_answer"
] is not True:
    raise SystemExit("the result dropped the provenance disclosure")
if payload["controls"][
    "lobular_necrosis_never_entered_a_composite_or_a_criterion"
] is not True:
    raise SystemExit("necrosis reached a composite")
if payload["controls"]["no_model_was_fitted"] is not True:
    raise SystemExit("a model was fitted in a gate that forbids one")

criteria = payload["criteria"]
d1 = criteria["d1_activity_and_fibrosis_are_independent"]["gate"]
d2 = criteria["d2_steatosis_and_saf_activity_are_independent"]["gate"]
if criteria["d1_activity_and_fibrosis_are_independent"]["decisive"] is not True:
    raise SystemExit("d1 lost its decisive mark")
for gate in (d1, d2):
    # met over applicable, never met over total, and never a silent all([]).
    if gate["n_conditions_applicable"] == 0 and gate["verdict"] != (
        "NO_APPLICABLE_CONDITIONS"
    ):
        raise SystemExit(f"{gate['gate']} passed with no applicable condition")
    if gate["passed"] != (
        gate["n_conditions_applicable"] > 0
        and gate["n_conditions_met"] == gate["n_conditions_applicable"]
    ):
        raise SystemExit(f"{gate['gate']} verdict does not follow from its conditions")
    if gate["n_conditions_applicable"] > gate["n_conditions_total"]:
        raise SystemExit(f"{gate['gate']} counts more applicable than total")

# the outcome must follow from the two checks
if d1["verdict"] == "NO_APPLICABLE_CONDITIONS" or d2["verdict"] == (
    "NO_APPLICABLE_CONDITIONS"
):
    expected = {"INDETERMINATE"}
elif d1["passed"] and d2["passed"]:
    expected = {"THREE_AXIS_SAF"}
elif d1["passed"]:
    expected = {"TWO_AXIS_ACTIVITY_FIBROSIS"}
elif d2["passed"]:
    expected = {"TWO_AXIS_STEATOSIS_ACTIVITY"}
else:
    expected = {"ONE_AXIS_ACTIVITY", "ONE_AXIS_FIBROSIS", "ONE_AXIS_UNDETERMINED"}
if payload["outcome"] not in expected:
    raise SystemExit("the outcome does not follow from the two gates")

for entry in payload["reuse_validation"]["association_agreement_with_masld_bench"]:
    if entry["max_absolute_difference_vs_masld_bench_spearman"] > 1e-10:
        raise SystemExit("the vectorized association disagrees with the package")
if payload["reuse_validation"]["bh_agreement_with_masld_bench"][
    "agrees_with_masld_bench_benjamini_hochberg"
] is not True:
    raise SystemExit("the BH step-up disagrees with the package")

# every criterion count must carry a measured floor
for name in criteria:
    for direction in criteria[name]["directions"]:
        if direction["count_null"]["n_draws"] < 1:
            raise SystemExit("a count was reported without a measured floor")

verify_frozen_tree(result)
freeze_tree(root, {
    "artifact_class": "gse267145_axis_count_campaign",
    "stage": "stage_0b_gate",
    "job": "2_of_2",
    "outcome": payload["outcome"],
    "prespecification_verified_before_computation": True,
    "parent_stage_0_decision": "STOP",
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "outcome": payload["outcome"],
    "grouping": payload["grouping"],
    "d1_verdict": d1["verdict"],
    "d1_met_over_applicable": f"{d1['n_conditions_met']}/{d1['n_conditions_applicable']}",
    "d2_verdict": d2["verdict"],
    "d2_met_over_applicable": f"{d2['n_conditions_met']}/{d2['n_conditions_applicable']}",
}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/axis_count/ARTIFACTS.json"
