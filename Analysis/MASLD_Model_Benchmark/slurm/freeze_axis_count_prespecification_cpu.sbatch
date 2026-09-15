#!/usr/bin/env bash
#SBATCH --job-name=model-check-1046
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 0b of the MASLD showcase model, job 1 of 2.
#
# Freezes what would count as one, two or three histological outcome axes
# BEFORE the counts that answer the question are computed. The analysis job
# (model-check-1047) verifies this tree's ARTIFACTS.json digest before it opens
# the expression matrix, so the ordering is provable from the output files alone
# and not from a claim in a report.
#
# The question was derived after seeing Stage 0's STOP. The prespecification
# says so in those words and cites Stage 0 by digest. What makes it a declared
# successor rather than a post hoc rationalisation is that its own criteria are
# fixed here, before its own data are opened.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/freeze_axis_count_prespecification.py
ENDPOINTS=${ROOT}/executions/model-data-061-21079623/activation/participant_endpoints.tsv
ENDPOINTS_SHA256=9a13c0fd65db0483566b7495def11ade3339861de091ea29fcbdc4c34f47b5b9
STAGE0_PRESPEC=${ROOT}/executions/model-check-1043-21181265-aspect-separability-prespec
STAGE0_RESULT=${ROOT}/executions/model-check-1044-21181334-aspect-separability
TARGET=${ROOT}/executions/model-check-1046-${SLURM_JOB_ID}-axis-count-prespec
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
bash -n slurm/freeze_axis_count_prespecification_cpu.sbatch

# The parent this question descends from must be the tree the prespecification
# cites, checked here rather than trusted from prose.
cat > "${STAGE}/stage0_inputs.sha256" <<'SUMS'
5395050a71a8e7a8bdc9ecd6b126c08e121a2f4ff6f42268bed0cbcdfe3efa76  executions/model-check-1043-21181265-aspect-separability-prespec/prespecification/ARTIFACTS.json
682238189f8bd97144eb7a7d39efefea8f17546d60d25ad655498af992f60c48  executions/model-check-1044-21181334-aspect-separability/separability/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/stage0_inputs.sha256"

# The endpoint table is the only other input this job reads. Its label geometry
# is what the determination audit is computed from, so the digest is checked
# before and after.
printf '%s  %s\n' "${ENDPOINTS_SHA256}" "${ENDPOINTS}" | sha256sum --check --strict

"${PYTHON}" -m unittest -v tests.unit.test_evaluate_gse267145_axis_count \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

# The script owns its --output path and refuses to overwrite one. The wrapper
# must not pre-create the directory it is handed.
"${PYTHON}" "${SCRIPT}" \
  --endpoints "${ENDPOINTS}" \
  --output "${STAGE}/prespecification/axis_count_prespecification.json" \
  > "${STAGE}/prespecification_stdout.json" 2> "${STAGE}/prespecification_stderr.txt"

printf '%s  %s\n' "${ENDPOINTS_SHA256}" "${ENDPOINTS}" | sha256sum --check --strict
sha256sum --check --strict "${STAGE}/stage0_inputs.sha256"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" \
  slurm/freeze_axis_count_prespecification_cpu.sbatch \
  tests/unit/test_evaluate_gse267145_axis_count.py \
  scripts/evaluate_gse267145_axis_count.py \
  "${ENDPOINTS}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
prespec_dir = root / "prespecification"
payload = json.loads(
    (prespec_dir / "axis_count_prespecification.json").read_text(encoding="utf-8")
)

expected = [
    "d1_activity_and_fibrosis_are_independent",
    "d2_steatosis_and_saf_activity_are_independent",
]
if sorted(payload["criteria"]) != expected:
    raise SystemExit("the two criteria are not the frozen two")
literal = (
    "in each direction, BH 0.05 count > 0 and above the 95th percentile of "
    "the matched within-sex permutation count null"
)
for name in expected:
    if payload["criteria"][name]["threshold"] != literal:
        raise SystemExit(f"{name} threshold differs from the frozen literal")
if payload["criteria"][expected[0]]["decisive"] is not True:
    raise SystemExit("d1 is not marked decisive")
if payload["criteria"][expected[1]]["decisive"] is not False:
    raise SystemExit("d2 must not be marked decisive")

provenance = payload["honest_provenance"]
if provenance["this_question_was_derived_after_seeing_stage_0s_answer"] is not True:
    raise SystemExit("the derivation after Stage 0 is not disclosed")
if provenance["criteria_fixed_before_any_expression_value_was_read"] is not True:
    raise SystemExit("the prespecification does not assert an unopened matrix")
if provenance["stage_0_prespecification_artifacts_sha256"] != (
    "5395050a71a8e7a8bdc9ecd6b126c08e121a2f4ff6f42268bed0cbcdfe3efa76"
):
    raise SystemExit("the cited Stage 0 digest is not the one checked above")
if payload["what_was_inspected_before_freezing"]["no_gene_expression_was_opened"] is not True:
    raise SystemExit("the prespecification does not assert an unopened matrix")

audit = payload["determination_audit_run_before_the_freeze"]
if audit["rank_fibrosis_in_the_span_of_the_three_aspects"] >= 0.60:
    raise SystemExit("fibrosis is inside the Stage 0 span; D1 would be a re-read")
if audit["stage0_did_not_deposit_a_fibrosis_or_composite_association_vector"] is not True:
    raise SystemExit("the audit does not record what Stage 0 deposited")

cells = payload["decision_rule"]["cells"]
if len(cells) != 4 or payload["decision_rule"]["outcome_map_is_complete_over_the_four_cells"] is not True:
    raise SystemExit("the outcome map does not cover the four cells")
if payload["diagnostics_never_gates"]["d3_does_a_composite_earn_its_place"]["is_a_gate"] is not False:
    raise SystemExit("d3 is not declared a diagnostic")

freeze_tree(prespec_dir, {
    "artifact_class": "gse267145_axis_count_prespecification",
    "prespec_id": payload["prespec_id"],
    "criteria_count": 2,
    "decisive_criterion": "d1_activity_and_fibrosis_are_independent",
    "derived_after_stage_0": True,
    "computed_before_any_axis_count": True,
    "gene_expression_opened": False,
    "status": "passed",
})
verify_frozen_tree(prespec_dir)
freeze_tree(root, {
    "artifact_class": "gse267145_axis_count_prespec_campaign",
    "stage": "stage_0b_gate",
    "job": "1_of_2",
    "parent_stage_0_decision": "STOP",
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "prespec_id": payload["prespec_id"],
    "criteria": sorted(payload["criteria"]),
    "outcomes": [cell["outcome"] for cell in cells],
    "determination_audit": {
        "fibrosis_in_the_aspect_span": audit[
            "rank_fibrosis_in_the_span_of_the_three_aspects"],
        "nas_sum_in_its_components_span": audit[
            "rank_nas_activity_sum_in_the_span_of_its_three_components"],
        "saf_activity_in_its_components_span": audit[
            "rank_saf_activity_sum_in_the_span_of_its_two_components"],
    },
}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/prespecification/ARTIFACTS.json"
