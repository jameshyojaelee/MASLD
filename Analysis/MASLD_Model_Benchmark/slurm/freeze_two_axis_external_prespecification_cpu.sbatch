#!/usr/bin/env bash
#SBATCH --job-name=model-check-1052
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 1, job 2 of 3: hold back what would count as external support.
#
# Held-back before any outcome-versus-expression association exists in either
# external arm. The substrate build (model-data-953) established the assignment
# and the join geometry and computed no external association, so this freeze is
# a real prespecification rather than a formality.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/freeze_two_axis_external_prespecification.py
SUBSTRATE_ROOT=${ROOT}/executions/model-data-953-21182514-two-axis-substrate
SUBSTRATE=${SUBSTRATE_ROOT}/substrate/substrate.json
ACTIVATION=${ROOT}/config/evaluation/two_axis_activity_fibrosis_external_activation.json
TARGET=${ROOT}/executions/model-check-1052-${SLURM_JOB_ID}-two-axis-external-prespec
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
bash -n slurm/freeze_two_axis_external_prespecification_cpu.sbatch

cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
68837a8293d5386e6bd7b539dd9107f61d1cfddb783feec373d9853e57be84ab  executions/model-data-953-21182514-two-axis-substrate/substrate/ARTIFACTS.json
5f26597ab6b6aeeaff867d54112554a33d4d0dc35be5f333e17e9010a7199fed  executions/model-check-1050-21182309-arm-b-axis-count/arm_b_axis_count/ARTIFACTS.json
58f72a4522309f9d4e0dfb1630c37fd994a308e714856b88508aa40cb27732f6  executions/model-check-1049-21182303-arm-b-axis-count-prespec/prespecification/ARTIFACTS.json
78ce955b6ef43ff45ebf78ec7744aa464febf3da93ce8aff4b6daf0f80c4ea3e  executions/model-check-1047-21181822-axis-count/axis_count/ARTIFACTS.json
91b18fbe8cb061f04825beb3d9c7651ade92e6055782eacd06b6ece9c2cab49f  executions/model-check-1048-21182175-axis-count-derived-after/addendum/ARTIFACTS.json
97db7a398b3b7dbb9a52c9e7599c3650b0c02b1277c5bcdc8431bde170fd267e  executions/model-check-1051-21182361-bootstrap-df-correction/correction/ARTIFACTS.json
3c5d3cc4e3f47b331c5f8ece04c213c304d49faa3ce3b974cef878c57e4f6ca4  config/evaluation/two_axis_activity_fibrosis_external_activation.json
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" - "${SUBSTRATE_ROOT}/substrate" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
print("substrate tree verified")
PY

"${PYTHON}" -m unittest -v tests.unit.test_freeze_two_axis_external_prespecification \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" "${SCRIPT}" \
  --substrate "${SUBSTRATE}" \
  --activation "${ACTIVATION}" \
  --output "${STAGE}/prespecification/two_axis_external_prespecification.json" \
  > "${STAGE}/prespecification_stdout.json" 2> "${STAGE}/prespecification_stderr.txt"

sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" \
  slurm/freeze_two_axis_external_prespecification_cpu.sbatch \
  tests/unit/test_freeze_two_axis_external_prespecification.py \
  "${SUBSTRATE}" "${ACTIVATION}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
prespec_dir = root / "prespecification"
payload = json.loads(
    (prespec_dir / "two_axis_external_prespecification.json").read_text(encoding="utf-8"))

expected = ["s1_the_assigned_set_is_elevated_on_its_own_axis",
            "s2_the_assigned_set_is_not_elevated_on_the_other_axis"]
if sorted(payload["criteria"]) != expected:
    raise SystemExit("the two criteria are not the frozen two")
if not payload["the_deliverable_and_its_honesty_constraint"][
        "external_support_is_set_level_not_per_gene"]:
    raise SystemExit("the set-level constraint is not declared")
unit = payload["unit_of_inference"]["GSE193066"]
if unit["n"] != 106 or unit["rows_in_the_deposited_matrix"] != 164:
    raise SystemExit("the GSE193066 unit is not 106 against a 164-row deposit")
if "NOT the analysis unit" not in unit["REQUIRED_NOTE_ON_THE_SUBSTRATE_ARTIFACT"]:
    raise SystemExit("the samples:164 note is missing")
if payload["unit_of_inference"]["GSE130970"]["distinct_people_assertable"] is not False:
    raise SystemExit("GSE130970 people were made assertable")
if payload["claim_boundary"][
    "the_fibrosis_axis_has_exactly_one_external_arm"] != (
        "The fibrosis axis has exactly one external arm, n=78 samples, "
        "with no assertable donor key."):
    raise SystemExit("the fibrosis claim boundary is not verbatim")
if payload["claim_boundary"]["activation_user_approved_scope_only"] is not False:
    raise SystemExit("the activation user_approved flag was not carried as recorded")
for name in expected:
    if payload["criteria"][name]["evaluated_on"] != ["activity_only", "fibrosis_only"]:
        raise SystemExit(f"{name} does not exclude the both cell")
if payload["diagnostics_never_gates"]["gse193066_fibrosis_arm"]["is_a_gate"] is not False:
    raise SystemExit("the GSE193066 fibrosis arm became a gate")
if payload["power"]["set_level_mde_is_the_operative_reference"]["is_a_gate"] is not False:
    raise SystemExit("the set-level MDE became a gate")
cells = [c["outcome"] for c in payload["decision_rule"]["cells"]]
for required in ("EXTERNALLY_SUPPORTED", "SUPPORTED_SINGLE_ARM",
                 "TRANSFERS_WITHOUT_SPECIFICITY", "NOT_EXTERNALLY_SUPPORTED",
                 "INDETERMINATE"):
    if required not in cells:
        raise SystemExit(f"the outcome map is missing {required}")

freeze_tree(prespec_dir, {
    "artifact_class": "two_axis_external_prespecification",
    "prespec_id": payload["prespec_id"],
    "criteria_count": 2,
    "external_support_is_set_level": True,
    "external_association_computed_before_the_freeze": False,
    "status": "passed",
})
verify_frozen_tree(prespec_dir)
freeze_tree(root, {
    "artifact_class": "two_axis_external_prespec_campaign",
    "stage": "stage_1", "job": "2_of_3", "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({"prespec_id": payload["prespec_id"], "outcomes": cells}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/prespecification/ARTIFACTS.json"
