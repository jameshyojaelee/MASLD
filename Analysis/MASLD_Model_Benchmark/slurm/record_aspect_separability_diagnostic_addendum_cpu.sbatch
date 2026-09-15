#!/usr/bin/env bash
#SBATCH --job-name=model-check-1045
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 0 of the MASLD showcase model, job 3 of 3.
#
# An additive overlay on the frozen prespecification (21181265) and the frozen
# result (21181334). It records the diagnostics' declared decision role and a
# post-hoc observation. It does NOT claim to have preceded the analysis: the
# instruction to freeze it first arrived after 21181334 had already completed,
# and the record says so. Both covered trees are verified by digest and must be
# byte-identical afterwards.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/record_aspect_separability_diagnostic_addendum.py

PRESPEC=${ROOT}/executions/model-check-1043-21181265-aspect-separability-prespec/prespecification
RESULT=${ROOT}/executions/model-check-1044-21181334-aspect-separability/separability
ENDPOINTS=${ROOT}/executions/model-data-061-21079623/activation/participant_endpoints.tsv
MOLECULAR=${ROOT}/executions/model-data-064-21079902/fixture/molecular

TARGET=${ROOT}/executions/model-check-1045-${SLURM_JOB_ID}-aspect-diagnostic-addendum
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
bash -n slurm/record_aspect_separability_diagnostic_addendum_cpu.sbatch

cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
5395050a71a8e7a8bdc9ecd6b126c08e121a2f4ff6f42268bed0cbcdfe3efa76  executions/model-check-1043-21181265-aspect-separability-prespec/prespecification/ARTIFACTS.json
31379d60091b3c77f32dfa8f43a1e75aefc702adf9cb7c874a115288f3319e8a  executions/model-check-1044-21181334-aspect-separability/ARTIFACTS.json
682238189f8bd97144eb7a7d39efefea8f17546d60d25ad655498af992f60c48  executions/model-check-1044-21181334-aspect-separability/separability/ARTIFACTS.json
4840ad08d9b2020fb6161856f4b3cccc054e38fce95b84653ccf7c1e440e14a2  executions/model-check-1044-21181334-aspect-separability/separability/aspect_separability.json
9a13c0fd65db0483566b7495def11ade3339861de091ea29fcbdc4c34f47b5b9  executions/model-data-061-21079623/activation/participant_endpoints.tsv
812d5fce6e7cc88f6939ae493db7b7e7470057c21cf3b219abc39f587ff64ddc  executions/model-data-064-21079902/fixture/molecular/rna_values.npy
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m unittest -v tests.unit.test_record_aspect_separability_diagnostic_addendum \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" "${SCRIPT}" \
  --root "${ROOT}" \
  --prespecification "${PRESPEC}" \
  --result "${RESULT}" \
  --molecular "${MOLECULAR}" \
  --endpoints "${ENDPOINTS}" \
  --analysis-job-id 21181334 \
  --analysis-end 2026-08-27T15:12:41 \
  --splits 200 \
  --seed 20260828 \
  --output "${STAGE}/addendum" \
  > "${STAGE}/addendum_stdout.json" 2> "${STAGE}/addendum_stderr.txt"

# The covered trees must be untouched by this overlay.
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" \
  slurm/record_aspect_separability_diagnostic_addendum_cpu.sbatch \
  tests/unit/test_record_aspect_separability_diagnostic_addendum.py \
  "${RESULT}/aspect_separability.json" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
record = root / "addendum"
addendum = json.loads(
    (record / "diagnostic_role_addendum.json").read_text(encoding="utf-8")
)
observations = json.loads(
    (record / "post_hoc_observations.json").read_text(encoding="utf-8")
)

ordering = addendum["honesty_about_ordering"]
if ordering["frozen_before_the_analysis"] is not False:
    raise SystemExit("the addendum must not claim it preceded the analysis")
if ordering["in_force_for_this_decision"] is not False:
    raise SystemExit("the addendum must not claim it governed this decision")
if addendum["frozen_verdict_altered"] is not False:
    raise SystemExit("the addendum altered a frozen verdict")
if addendum["promotion_gate_registered"] is not False:
    raise SystemExit("the addendum registered a gate")
role = addendum["declared_decision_role"]
if role["may_downgrade_a_go_to_indeterminate"] is not True:
    raise SystemExit("the declared role lost its downgrade direction")
if role["may_never_upgrade_a_stop_to_a_go"] is not True:
    raise SystemExit("the declared role lost its asymmetry")
if observations["explicitly_post_hoc"] is not True:
    raise SystemExit("the observation is not marked post hoc")
if observations["authorises_nothing"] is not True:
    raise SystemExit("the observation claims authority")

verify_frozen_tree(record)
freeze_tree(root, {
    "artifact_class": "gse267145_aspect_diagnostic_addendum_campaign",
    "stage": "stage_0_gate",
    "job": "3_of_3",
    "record_is_additive_overlay": True,
    "frozen_verdict_altered": False,
    "in_force_for_the_frozen_decision": False,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "frozen_before_the_analysis": ordering["frozen_before_the_analysis"],
    "in_force_for_this_decision": ordering["in_force_for_this_decision"],
    "could_it_have_changed_this_decision": ordering[
        "could_it_have_changed_this_decision"],
    "proposed_ceiling_would_reject_fibrosis": addendum["reliability_gate"][
        "measured_counterexample"]["observed_exceeds_the_proposed_ceiling"],
}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/addendum/ARTIFACTS.json"
