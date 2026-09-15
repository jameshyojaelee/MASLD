#!/usr/bin/env bash
#SBATCH --job-name=model-check-1054
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 3, job 2 of 3: hold back what would count as composition-independent
# association. Held-back before any composition-adjusted number exists.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/freeze_composition_prespecification.py
SUB_ROOT=${ROOT}/executions/model-data-954-21183077-composition-substrate/substrate
SUBSTRATE=${SUB_ROOT}/composition_substrate.json
TARGET=${ROOT}/executions/model-check-1054-${SLURM_JOB_ID}-composition-prespec
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
bash -n slurm/freeze_composition_prespecification_cpu.sbatch

cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
97f6c12407ac5d093afc978a777983b14804d971963c3291d6510059f9112623  executions/model-data-954-21183077-composition-substrate/substrate/ARTIFACTS.json
daf4f02db0585f96b5605c09f4f71d6182e77d2d9bae25e263b34b1ba2e7eece  executions/model-check-1053-21182545-two-axis-external/external/ARTIFACTS.json
5f26597ab6b6aeeaff867d54112554a33d4d0dc35be5f333e17e9010a7199fed  executions/model-check-1050-21182309-arm-b-axis-count/arm_b_axis_count/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

"${PYTHON}" - "${SUB_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]).resolve(strict=True))
print("composition substrate verified")
PY

"${PYTHON}" -m unittest -v tests.unit.test_freeze_composition_prespecification \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" "${SCRIPT}" --substrate "${SUBSTRATE}" \
  --output "${STAGE}/prespecification/composition_prespecification.json" \
  > "${STAGE}/prespecification_stdout.json" 2> "${STAGE}/prespecification_stderr.txt"

sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" slurm/freeze_composition_prespecification_cpu.sbatch \
  tests/unit/test_freeze_composition_prespecification.py \
  "${SUBSTRATE}" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
d = root / "prespecification"
p = json.loads((d / "composition_prespecification.json").read_text(encoding="utf-8"))
if sorted(p["criteria"]) != ["t1_a_composition_independent_component_exists",
                             "t2_per_gene_composition_dependence"]:
    raise SystemExit("the two criteria are not the frozen two")
if p["criteria"]["t1_a_composition_independent_component_exists"][
        "failure_is_named"] != "SET_IS_COMPOSITION_MEDIATED":
    raise SystemExit("T1's failure is not named for mechanism")
if not p["the_framing_that_governs_every_reading"][
        "a_composition_dependent_gene_is_not_an_artifact"]:
    raise SystemExit("the mediation framing is missing")
c = p["correction_to_the_frozen_substrate_artifact"]
if "three are stable, not two" not in c["what_is_true"]:
    raise SystemExit("the variant-stability correction is missing or wrong")
if len(p["substrate"]["eligible_family"]) != 6:
    raise SystemExit("the family is not the approved six")
if p["covariates"]["n_covariates"] != 7:
    raise SystemExit("the covariate count is not seven")
for arm in p["covariates"]["conditioning_is_a_measured_negative"][
        "measured_per_arm"].values():
    if arm["numerical_rank"] != arm["n_covariates"]:
        raise SystemExit("a covariate block is rank deficient")
cells = [x["outcome"] for x in p["decision_rule"]["cells"]]
if "SET_IS_COMPOSITION_MEDIATED" not in cells or "INDETERMINATE" not in cells:
    raise SystemExit("the outcome map is incomplete")
freeze_tree(d, {"artifact_class": "composition_prespecification",
                "prespec_id": p["prespec_id"], "criteria_count": 2,
                "family_size": 6, "composition_adjusted_number_existed": False,
                "status": "passed"})
verify_frozen_tree(d)
freeze_tree(root, {"artifact_class": "composition_prespec_campaign",
                   "stage": "stage_3", "job": "2_of_3", "status": "passed"})
verify_frozen_tree(root)
print(json.dumps({"prespec_id": p["prespec_id"], "outcomes": cells}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/prespecification/ARTIFACTS.json"
