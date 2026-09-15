#!/usr/bin/env bash
#SBATCH --job-name=model-check-1055
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 3, job 3 of 3: the composition-dependence check.
# The prespecification frozen by model-check-1054 (21183086) is verified before
# any composition-adjusted association is computed.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/evaluate_composition_dependence.py
SEALED=${ROOT}/scripts/evaluate_gse267145_axis_count.py
EXT=${ROOT}/scripts/multicovariate_partial.py
PRESPEC_ROOT=${ROOT}/executions/model-check-1054-21183086-composition-prespec
PRESPEC=${PRESPEC_ROOT}/prespecification
STAGE1=${ROOT}/executions/model-data-953-21182514-two-axis-substrate/substrate
COMP=${ROOT}/executions/model-data-954-21183077-composition-substrate/substrate
ARM_B=${ROOT}/executions/model-data-880-21130257-gse135251-source/source
EXTERNAL=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk

TARGET=${ROOT}/executions/model-check-1055-${SLURM_JOB_ID}-composition-dependence
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the composition result.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/evaluate_composition_dependence_cpu.sbatch

cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
4e837504d7a7a4968b99f04276a4c232a3c5b84a2c9a6c253725a3f9fe1f83f5  executions/model-check-1054-21183086-composition-prespec/prespecification/ARTIFACTS.json
97f6c12407ac5d093afc978a777983b14804d971963c3291d6510059f9112623  executions/model-data-954-21183077-composition-substrate/substrate/ARTIFACTS.json
68837a8293d5386e6bd7b539dd9107f61d1cfddb783feec373d9853e57be84ab  executions/model-data-953-21182514-two-axis-substrate/substrate/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"
printf '%s  %s\n' "5f4e483dc22f1143f1c4447f60e308a5822dd86070540527abbd429e7278e559" "${SEALED}" | sha256sum --check --strict

"${PYTHON}" - "${PRESPEC}" "${PRESPEC_ROOT}" "${COMP}" "${STAGE1}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for a in sys.argv[1:]:
    verify_frozen_tree(Path(a).resolve(strict=True))
print("frozen trees verified")
PY

"${PYTHON}" -m unittest -v tests.unit.test_evaluate_composition_dependence \
  tests.unit.test_multicovariate_partial > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" "${SCRIPT}" \
  --prespecification "${PRESPEC}" --sealed-instrument "${SEALED}" \
  --extension "${EXT}" --stage1-substrate "${STAGE1}" \
  --composition-substrate "${COMP}" --arm-b-source "${ARM_B}" \
  --external-root "${EXTERNAL}" --output "${STAGE}/composition" \
  --null-draws 2000 --mde-draws 400 --seed 20260827 \
  > "${STAGE}/composition_stdout.txt" 2> "${STAGE}/composition_stderr.txt"

sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" "${SEALED}" "${EXT}" \
  slurm/evaluate_composition_dependence_cpu.sbatch \
  tests/unit/test_evaluate_composition_dependence.py \
  "${PRESPEC}/composition_prespecification.json" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
r = root / "composition"
p = json.loads((r / "composition_dependence.json").read_text(encoding="utf-8"))
declared = {"COMPOSITION_INDEPENDENT_COMPONENT_EXISTS",
            "COMPOSITION_INDEPENDENT_IN_ONE_ARM_ONLY",
            "SET_IS_COMPOSITION_MEDIATED", "INDETERMINATE"}
if p["outcome"] not in declared:
    raise SystemExit(f"undeclared outcome {p['outcome']!r}")
if not p["the_framing_that_governs_every_reading"][
        "a_composition_dependent_gene_is_not_an_artifact"]:
    raise SystemExit("the mediation framing was dropped")
if len(p["eligible_family"]) != 6:
    raise SystemExit("the family is not the approved six")
g = p["criteria"]["t1_a_composition_independent_component_exists"]["gate"]
if g["n_conditions_applicable"] == 0 and g["verdict"] != "NO_APPLICABLE_CONDITIONS":
    raise SystemExit("T1 passed with no applicable condition")
if g["passed"] != (g["n_conditions_applicable"] > 0
                   and g["n_conditions_met"] == g["n_conditions_applicable"]):
    raise SystemExit("the T1 verdict does not follow from its conditions")
for c in g["conditions"]:
    if c.get("verdict") == "set_is_composition_mediated" and not c["applicable"]:
        raise SystemExit("a mediated verdict was recorded on an inapplicable condition")
for cohort, arm in p["arms"].items():
    for axis, cond in arm["conditioning"].items():
        if cond["numerical_rank"] != cond["n_covariates"]:
            raise SystemExit(f"{cohort}/{axis} covariate block is rank deficient")
if p["arms"]["GSE193066"]["n"] != 106:
    raise SystemExit("GSE193066 was not restricted to 106 first-biopsy participants")
if p["arms"]["GSE193066"]["axes_served"] != ["activity"]:
    raise SystemExit("GSE193066 served an axis it is not scoped for")
verify_frozen_tree(r)
freeze_tree(root, {"artifact_class": "composition_dependence_campaign",
                   "stage": "stage_3", "job": "3_of_3", "outcome": p["outcome"],
                   "prespecification_verified_before_computation": True,
                   "status": "passed"})
verify_frozen_tree(root)
print(json.dumps({"outcome": p["outcome"]}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/composition/ARTIFACTS.json"
