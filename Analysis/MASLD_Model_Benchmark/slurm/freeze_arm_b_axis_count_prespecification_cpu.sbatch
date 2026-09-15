#!/usr/bin/env bash
#SBATCH --job-name=model-check-1049
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 0c of the MASLD showcase model, job 1 of 2.
#
# Holds back whether activity and fibrosis are two axes in GSE135251 at n=180,
# BEFORE any number has been computed in that cohort. That restraint is the
# whole point: Arm B is untouched, so holding back first buys a genuine
# prespecification rather than a formality. The freeze reads only the two axes'
# marginals and tie structure and never their correlation.
#
# The parents this question descends from -- Stage 0's STOP and Stage 0b's
# ONE_AXIS_ACTIVITY -- are pinned by digest and their numbers are recorded by
# value inside the frozen document, so the derived-after exposure is on the
# record rather than left for a reader to reconstruct.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/freeze_arm_b_axis_count_prespecification.py
ARMB=${ROOT}/executions/model-data-880-21130257-gse135251-source/source
ENDPOINTS=${ARMB}/outcomes/participant_endpoints.tsv
TARGET=${ROOT}/executions/model-check-1049-${SLURM_JOB_ID}-arm-b-axis-count-prespec
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the prespecification record.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/freeze_arm_b_axis_count_prespecification_cpu.sbatch

# Every parent, and the untouched Arm B substrate, checked rather than trusted.
cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
5395050a71a8e7a8bdc9ecd6b126c08e121a2f4ff6f42268bed0cbcdfe3efa76  executions/model-check-1043-21181265-aspect-separability-prespec/prespecification/ARTIFACTS.json
4840ad08d9b2020fb6161856f4b3cccc054e38fce95b84653ccf7c1e440e14a2  executions/model-check-1044-21181334-aspect-separability/separability/aspect_separability.json
8775349ea96abfa9563318f873f360dd6cbe2391ca34f571f9a3100d2851e024  executions/model-check-1045-21181430-aspect-diagnostic-addendum/addendum/post_hoc_observations.json
e7397c3038dfbbdd1619048acfdd6cabc0e018a7af161175925130a7a6ffed3e  executions/model-check-1046-21181811-axis-count-prespec/prespecification/ARTIFACTS.json
78ce955b6ef43ff45ebf78ec7744aa464febf3da93ce8aff4b6daf0f80c4ea3e  executions/model-check-1047-21181822-axis-count/axis_count/ARTIFACTS.json
91b18fbe8cb061f04825beb3d9c7651ade92e6055782eacd06b6ece9c2cab49f  executions/model-check-1048-21182175-axis-count-derived-after/addendum/ARTIFACTS.json
9e7bdda3d4511a9bc0e520515f9ab12f4ca0729ca9a74a33f2c0491a4f253db7  executions/model-data-880-21130257-gse135251-source/source/outcomes/participant_endpoints.tsv
d6219f67fa676e888f47b0225a75a18036ddf83ba9bbd2f1991069e3781e39aa  executions/model-data-880-21130257-gse135251-source/source/molecular/rna_values.npy
3dcdd51055010215592f9ff7372c62824b7345565d1bda4428a29536b5e5c3a8  executions/model-data-880-21130257-gse135251-source/source/molecular/rna_feature_axis.tsv
8a959b3da6feafe419170559c9f3ac92803c3d683c8818a60a4da09ab33f619b  executions/model-data-880-21130257-gse135251-source/source/molecular/participant_axis.tsv
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m unittest -v tests.unit.test_freeze_arm_b_axis_count_prespecification \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

# The expression matrix must be untouched by this job. Its digest is checked
# above and again below, and the freeze script is never handed its path.
"${PYTHON}" "${SCRIPT}" \
  --endpoints "${ENDPOINTS}" \
  --output "${STAGE}/prespecification/arm_b_axis_count_prespecification.json" \
  > "${STAGE}/prespecification_stdout.json" 2> "${STAGE}/prespecification_stderr.txt"

sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" \
  slurm/freeze_arm_b_axis_count_prespecification_cpu.sbatch \
  tests/unit/test_freeze_arm_b_axis_count_prespecification.py \
  "${ENDPOINTS}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
prespec_dir = root / "prespecification"
payload = json.loads(
    (prespec_dir / "arm_b_axis_count_prespecification.json").read_text(encoding="utf-8")
)

expected = [
    "e1_activity_and_fibrosis_are_independent",
    "e2_an_empty_direction_is_classified_against_the_measured_floor",
]
if sorted(payload["criteria"]) != expected:
    raise SystemExit("the two criteria are not the frozen two")
if payload["criteria"][expected[0]]["decisive"] is not True:
    raise SystemExit("e1 is not marked decisive")
if payload["criteria"][expected[1]]["decisive"] is not False:
    raise SystemExit("e2 must not be marked decisive")
if payload["criteria"][expected[1]][
    "is_a_disambiguation_rule_not_a_pass_fail_gate"] is not True:
    raise SystemExit("e2 is not marked a disambiguation rule")

provenance = payload["honest_provenance"]
if provenance["this_question_was_derived_after_two_completed_results"] is not True:
    raise SystemExit("the derivation after Stage 0 and 0b is not disclosed")
if provenance["criteria_fixed_before_any_arm_b_number_was_computed"] is not True:
    raise SystemExit("the prespecification does not assert an untouched Arm B")
if "in-cohort" not in provenance["in_cohort_and_why_that_is_acceptable"]:
    raise SystemExit("the in-cohort status is not stated")

inspected = payload["what_was_inspected_before_freezing"]
if inspected["no_gene_expression_was_opened"] is not True:
    raise SystemExit("the prespecification does not assert an unopened matrix")
if inspected["no_correlation_between_the_two_axes_was_computed"] is not True:
    raise SystemExit("the restraint on the cross-axis correlation is not asserted")

# The restraint must hold in the document, not only be claimed by it.
blob = json.dumps(payload)
for forbidden in ("nas_score|fibrosis_stage", "fibrosis_stage|nas_score"):
    if forbidden in blob:
        raise SystemExit(f"the document carries a cross-axis number: {forbidden}")

hazard = payload["fibrosis_scale_hazard"]
if "not determinable" not in hazard[
    "whether_arm_a_lacks_a_level_4_or_sampled_no_f4_participant"]:
    raise SystemExit("the scale question was guessed rather than left undetermined")
if payload["controls"]["no_cross_substrate_fibrosis_comparison_is_made"] is not True:
    raise SystemExit("a cross-substrate fibrosis comparison was declared")

cells = [cell["outcome"] for cell in payload["decision_rule"]["cells"]]
for required in (
    "TWO_AXIS_ACTIVITY_FIBROSIS",
    "ONE_AXIS_ACTIVITY_FIBROSIS_UNDERPOWERED",
    "ONE_AXIS_ACTIVITY_FIBROSIS_NOT_INDEPENDENT",
    "ONE_AXIS_FIBROSIS_ACTIVITY_UNDERPOWERED",
    "ONE_AXIS_FIBROSIS_ACTIVITY_NOT_INDEPENDENT",
    "NO_AXIS_RESOLVED",
):
    if required not in cells:
        raise SystemExit(f"the outcome map is missing {required}")

prior = payload["arm_a_prior_observations"]
if prior["partial_families_bh_0_05_counts"]["nas_activity_sum|fibrosis"] != 805:
    raise SystemExit("the Arm A observations are not recorded by value")
if payload["cohort"]["n_participants"] != 180:
    raise SystemExit("the participant count is not 180")

freeze_tree(prespec_dir, {
    "artifact_class": "gse135251_axis_count_prespecification",
    "prespec_id": payload["prespec_id"],
    "criteria_count": 2,
    "decisive_criterion": "e1_activity_and_fibrosis_are_independent",
    "derived_after_stage_0_and_stage_0b": True,
    "arm_b_untouched_before_the_freeze": True,
    "gene_expression_opened": False,
    "cross_axis_correlation_computed": False,
    "status": "passed",
})
verify_frozen_tree(prespec_dir)
freeze_tree(root, {
    "artifact_class": "gse135251_axis_count_prespec_campaign",
    "stage": "stage_0c_gate",
    "job": "1_of_2",
    "parent_stage_0b_outcome": "ONE_AXIS_ACTIVITY",
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "prespec_id": payload["prespec_id"],
    "n": payload["cohort"]["n_participants"],
    "criteria": sorted(payload["criteria"]),
    "outcomes": cells,
}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/prespecification/ARTIFACTS.json"
