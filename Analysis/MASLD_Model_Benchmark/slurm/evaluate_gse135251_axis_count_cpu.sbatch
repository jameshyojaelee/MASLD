#!/usr/bin/env bash
#SBATCH --job-name=model-check-1050
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 0c, job 2 of 2: the Arm B axis-count check at n=180.
#
# The prespecification frozen by model-check-1049 (job 21182303) is verified by
# digest BEFORE any GSE135251 number is computed. Arm B was untouched at the
# freeze, so this ordering is what makes the prespecification real rather than a
# formality. Every statistic is imported from the Stage 0b analysis script, so
# the two substrates are compared on one instrument.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/evaluate_gse135251_axis_count.py
ARM_A=${ROOT}/scripts/evaluate_gse267145_axis_count.py
ADDENDUM=${ROOT}/scripts/record_axis_count_derived_after_addendum.py
PRESPEC_ROOT=${ROOT}/executions/model-check-1049-21182303-arm-b-axis-count-prespec
PRESPEC=${PRESPEC_ROOT}/prespecification
ARMB=${ROOT}/executions/model-data-880-21130257-gse135251-source/source
ENDPOINTS=${ARMB}/outcomes/participant_endpoints.tsv
MOLECULAR=${ARMB}/molecular

TARGET=${ROOT}/executions/model-check-1050-${SLURM_JOB_ID}-arm-b-axis-count
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the Arm B axis-count result.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/evaluate_gse135251_axis_count_cpu.sbatch

cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
181793b9e3d9d17adfe590c4ae6ed73f7f495a4eeb3f54b324611f400d241f1c  executions/model-check-1049-21182303-arm-b-axis-count-prespec/ARTIFACTS.json
58f72a4522309f9d4e0dfb1630c37fd994a308e714856b88508aa40cb27732f6  executions/model-check-1049-21182303-arm-b-axis-count-prespec/prespecification/ARTIFACTS.json
78ce955b6ef43ff45ebf78ec7744aa464febf3da93ce8aff4b6daf0f80c4ea3e  executions/model-check-1047-21181822-axis-count/axis_count/ARTIFACTS.json
91b18fbe8cb061f04825beb3d9c7651ade92e6055782eacd06b6ece9c2cab49f  executions/model-check-1048-21182175-axis-count-derived-after/addendum/ARTIFACTS.json
9e7bdda3d4511a9bc0e520515f9ab12f4ca0729ca9a74a33f2c0491a4f253db7  executions/model-data-880-21130257-gse135251-source/source/outcomes/participant_endpoints.tsv
d6219f67fa676e888f47b0225a75a18036ddf83ba9bbd2f1991069e3781e39aa  executions/model-data-880-21130257-gse135251-source/source/molecular/rna_values.npy
3dcdd51055010215592f9ff7372c62824b7345565d1bda4428a29536b5e5c3a8  executions/model-data-880-21130257-gse135251-source/source/molecular/rna_feature_axis.tsv
8a959b3da6feafe419170559c9f3ac92803c3d683c8818a60a4da09ab33f619b  executions/model-data-880-21130257-gse135251-source/source/molecular/participant_axis.tsv
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" - "${PRESPEC}" "${PRESPEC_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for argument in sys.argv[1:]:
    verify_frozen_tree(Path(argument).resolve(strict=True))
print("prespecification tree verified")
PY

"${PYTHON}" -m unittest -v tests.unit.test_evaluate_gse135251_axis_count \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

find "${MOLECULAR}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_before.txt"

"${PYTHON}" "${SCRIPT}" \
  --prespecification "${PRESPEC}" \
  --analysis-script "${ARM_A}" \
  --addendum-script "${ADDENDUM}" \
  --endpoints "${ENDPOINTS}" \
  --molecular "${MOLECULAR}" \
  --output "${STAGE}/arm_b_axis_count" \
  --count-null-draws 20000 \
  --null-block 250 \
  --splits 200 \
  --seed 20260827 \
  > "${STAGE}/arm_b_stdout.txt" 2> "${STAGE}/arm_b_stderr.txt"

find "${MOLECULAR}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_after.txt"
cmp -s "${STAGE}/frozen_before.txt" "${STAGE}/frozen_after.txt" || { printf 'A frozen artifact changed.\n' >&2; exit 1; }
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" "${ARM_A}" "${ADDENDUM}" \
  slurm/evaluate_gse135251_axis_count_cpu.sbatch \
  tests/unit/test_evaluate_gse135251_axis_count.py \
  "${PRESPEC}/arm_b_axis_count_prespecification.json" \
  "${ENDPOINTS}" "${MOLECULAR}/rna_values.npy" \
  "${MOLECULAR}/rna_feature_axis.tsv" \
  "${MOLECULAR}/participant_axis.tsv" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
result = root / "arm_b_axis_count"
payload = json.loads((result / "arm_b_axis_count.json").read_text(encoding="utf-8"))

if payload["prespec_id"] != "gse135251_axis_count_v1":
    raise SystemExit("the result is not against the frozen prespecification")
declared = {
    "TWO_AXIS_ACTIVITY_FIBROSIS",
    "ONE_AXIS_ACTIVITY_FIBROSIS_UNDERPOWERED",
    "ONE_AXIS_ACTIVITY_FIBROSIS_NOT_INDEPENDENT",
    "ONE_AXIS_FIBROSIS_ACTIVITY_UNDERPOWERED",
    "ONE_AXIS_FIBROSIS_ACTIVITY_NOT_INDEPENDENT",
    "NO_AXIS_RESOLVED",
    "INDETERMINATE",
}
if payload["outcome"] not in declared:
    raise SystemExit(f"undeclared outcome {payload['outcome']!r}")
if payload["cohort"]["distinct_participants"] != 180:
    raise SystemExit("the participant count is not 180")
if payload["controls"]["no_model_is_fitted"] is not True:
    raise SystemExit("a model was fitted in a gate that forbids one")
if payload["controls"]["no_cross_substrate_fibrosis_comparison_is_made"] is not True:
    raise SystemExit("a cross-substrate fibrosis comparison was made")

criteria = payload["criteria"]
e1 = criteria["e1_activity_and_fibrosis_are_independent"]["gate"]
e2 = criteria["e2_an_empty_direction_is_classified_against_the_measured_floor"]
if criteria["e1_activity_and_fibrosis_are_independent"]["decisive"] is not True:
    raise SystemExit("e1 lost its decisive mark")
if e2["decisive"] is not False:
    raise SystemExit("e2 became decisive")
# met over applicable, never met over total, and never a silent all([]).
if e1["n_conditions_applicable"] == 0 and e1["verdict"] != "NO_APPLICABLE_CONDITIONS":
    raise SystemExit("e1 passed with no applicable condition")
if e1["passed"] != (
    e1["n_conditions_applicable"] > 0
    and e1["n_conditions_met"] == e1["n_conditions_applicable"]
):
    raise SystemExit("the e1 verdict does not follow from its conditions")

# the outcome must follow from E1, qualified by E2 and never overridden by it
if e1["verdict"] == "NO_APPLICABLE_CONDITIONS":
    expected = {"INDETERMINATE"}
elif e1["passed"]:
    expected = {"TWO_AXIS_ACTIVITY_FIBROSIS"}
else:
    survivors = [c for c in e1["conditions"] if c["applicable"] and c["met"]]
    if not survivors:
        expected = {"NO_AXIS_RESOLVED"}
    else:
        expected = {
            "ONE_AXIS_ACTIVITY_FIBROSIS_UNDERPOWERED",
            "ONE_AXIS_ACTIVITY_FIBROSIS_NOT_INDEPENDENT",
            "ONE_AXIS_FIBROSIS_ACTIVITY_UNDERPOWERED",
            "ONE_AXIS_FIBROSIS_ACTIVITY_NOT_INDEPENDENT",
            "INDETERMINATE",
        }
if payload["outcome"] not in expected:
    raise SystemExit("the outcome does not follow from the gates")

for entry in payload["reuse_validation"]["association_agreement_with_masld_bench"]:
    if entry["max_absolute_difference_vs_masld_bench_spearman"] > 1e-10:
        raise SystemExit("the vectorized association disagrees with the package")
for direction in criteria["e1_activity_and_fibrosis_are_independent"]["directions"]:
    if direction["count_null"]["n_draws"] < 1:
        raise SystemExit("a count was reported without a measured floor")

verify_frozen_tree(result)
freeze_tree(root, {
    "artifact_class": "gse135251_axis_count_campaign",
    "stage": "stage_0c_gate",
    "job": "2_of_2",
    "outcome": payload["outcome"],
    "prespecification_verified_before_computation": True,
    "arm_b_untouched_before_the_freeze": True,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "outcome": payload["outcome"],
    "e1_verdict": e1["verdict"],
    "e1_met_over_applicable":
        f"{e1['n_conditions_met']}/{e1['n_conditions_applicable']}",
}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/arm_b_axis_count/ARTIFACTS.json"
