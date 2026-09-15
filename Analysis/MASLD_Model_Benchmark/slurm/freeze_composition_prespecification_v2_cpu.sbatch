#!/usr/bin/env bash
#SBATCH --job-name=model-check-1060
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 3 v2, job 1 of 2: hold back again the composition criteria after three defects.
#
# Written by hand. Deriving this wrapper from the v1 one by sed produced three
# separate provenance defects across two attempts -- a doubled _v2 suffix, a
# dotted module path turned into a filename, and a manifest hashing the v1
# sibling that did not run. The manifest and unit-test lines are the ones that
# carry the provenance, so they are typed rather than transformed.
#
# This supersedes the hold-back from job 21183286, whose source manifest named the v1
# sbatch and whose unit-test check ran the evaluator's tests rather than any test
# of the freeze itself.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/freeze_composition_prespecification_v2.py
WRAPPER=slurm/freeze_composition_prespecification_v2_cpu.sbatch
FREEZER_TEST=tests/unit/test_freeze_composition_prespecification_v2.py
SUB_ROOT=${ROOT}/executions/model-data-954-21183077-composition-substrate/substrate
SUBSTRATE=${SUB_ROOT}/composition_substrate.json

TARGET=${ROOT}/executions/model-check-1060-${SLURM_JOB_ID}-composition-prespec-v2
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the prespecification record.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n "${WRAPPER}"

# Every file the manifest will name must exist before the job does any work,
# so a mistyped path fails here rather than after the output file is written.
for required in "${SCRIPT}" "${WRAPPER}" "${FREEZER_TEST}" "${SUBSTRATE}"; do
  [[ -f "${required}" ]] || { printf 'Missing: %s\n' "${required}" >&2; exit 1; }
done

cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
97f6c12407ac5d093afc978a777983b14804d971963c3291d6510059f9112623  executions/model-data-954-21183077-composition-substrate/substrate/ARTIFACTS.json
4e837504d7a7a4968b99f04276a4c232a3c5b84a2c9a6c253725a3f9fe1f83f5  executions/model-check-1054-21183086-composition-prespec/prespecification/ARTIFACTS.json
e6e9bb04064fe1f1d9df1f8346b6633f819e6953ccb00d2e4c04ba6886d610cb  executions/model-check-1055-21183164-composition-dependence/composition/ARTIFACTS.json
68837a8293d5386e6bd7b539dd9107f61d1cfddb783feec373d9853e57be84ab  executions/model-data-953-21182514-two-axis-substrate/substrate/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" - "${SUB_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
print("composition substrate verified")
PY

# The freeze is conditional by a test OF THE FREEZE. The v1 v2-attempt ran the
# evaluator's tests here, which pinned nothing about the prespecification.
"${PYTHON}" -m unittest -v tests.unit.test_freeze_composition_prespecification_v2 \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" "${SCRIPT}" --substrate "${SUBSTRATE}" \
  --output "${STAGE}/prespecification/composition_prespecification_v2.json" \
  > "${STAGE}/prespecification_stdout.json" 2> "${STAGE}/prespecification_stderr.txt"

sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"

# The manifest records what actually ran: this wrapper, this script, this test.
sha256sum "${SCRIPT}" "${WRAPPER}" "${FREEZER_TEST}" "${SUBSTRATE}" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${WRAPPER}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
wrapper = sys.argv[2]
manifest = (root / "source.sha256").read_text(encoding="utf-8")
if wrapper not in manifest:
    raise SystemExit(f"the source manifest does not name the wrapper that ran: {wrapper}")
if "freeze_composition_prespecification_cpu.sbatch" in manifest:
    raise SystemExit("the source manifest still names the v1 sbatch, which did not run")
if "test_freeze_composition_prespecification_v2.py" not in manifest:
    raise SystemExit("the source manifest does not name the freezer's own test")

d = root / "prespecification"
p = json.loads((d / "composition_prespecification_v2.json").read_text(encoding="utf-8"))
if p["prespec_id"] != "composition_dependence_v2":
    raise SystemExit("this is not the v2 prespecification")
if p["supersedes"] != "composition_dependence_v1":
    raise SystemExit("the superseded prespecification is not recorded")
prior = p["the_prior_run_is_not_a_result"]
if prior["may_it_be_cited"] is not False:
    raise SystemExit("the void run was not marked uncitable")
if prior["cause_1_the_training_arm_dropped_out_silently"]["measured"][
        "join_on_feature_index"] != 0:
    raise SystemExit("the measured cause of the silent dropout is missing")
t1 = p["criteria"]["t1_a_composition_independent_component_exists"]
if "BY AT LEAST" not in t1["threshold"] or "either side" not in t1["threshold"]:
    raise SystemExit("T1's guard is not symmetric in the frozen literal")
if t1["failure_is_named"] != "SET_IS_COMPOSITION_MEDIATED":
    raise SystemExit("T1's failure name drifted")
t2 = p["criteria"]["t2_per_gene_composition_dependence"]
if t2["confirmatory_arm"] != "GSE135251":
    raise SystemExit("T2 is not scoped to the confirmatory arm")
if set(t2["external_arms_are_indeterminate_by_construction"]["arms"]) != {
        "GSE130970", "GSE193066"}:
    raise SystemExit("the external arms are not scoped indeterminate")
cells = [x["outcome"] for x in p["decision_rule"]["cells"]]
if "ARMS_DIFFER_WITHIN_RESOLUTION" not in cells:
    raise SystemExit("the within-resolution cell is missing")
if len(p["substrate"]["eligible_family"]) != 6:
    raise SystemExit("the family is not the approved six")
if "never trust the position" not in p[
        "the_column_lesson_runs_both_ways_in_this_tree"]["from_this_defect"]:
    raise SystemExit("the column lesson is not recorded in both directions")

freeze_tree(d, {"artifact_class": "composition_prespecification_v2",
                "prespec_id": p["prespec_id"], "criteria_count": 2,
                "supersedes": "composition_dependence_v1",
                "family_size": 6, "status": "passed"})
verify_frozen_tree(d)
freeze_tree(root, {"artifact_class": "composition_prespec_v2_campaign",
                   "stage": "stage_3_v2", "job": "1_of_2",
                   "manifest_names_the_wrapper_that_ran": True,
                   "gated_by_a_test_of_the_freeze": True,
                   "status": "passed"})
verify_frozen_tree(root)
print(json.dumps({"prespec_id": p["prespec_id"], "outcomes": cells}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/prespecification/ARTIFACTS.json"
