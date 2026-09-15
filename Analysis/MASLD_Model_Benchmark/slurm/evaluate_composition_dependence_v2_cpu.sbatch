#!/usr/bin/env bash
#SBATCH --job-name=model-check-1061
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 3 v2, job 2 of 2: the composition-dependence check.
#
# Written by hand, like its freeze counterpart. Deriving a wrapper from a
# sibling by sed produced three provenance defects on this stage; the manifest
# and unit-test lines carry the provenance and are typed, not transformed.
#
# Binds the canonical v2 prespecification, model-check-1058-21183309. A second
# tree, 1057-21183286, holds a byte-identical specification body but a source
# manifest naming a wrapper that did not run; the analysis records that as a
# duplicate hold back rather than leaving two objects claiming one identity.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/evaluate_composition_dependence_v2.py
WRAPPER=slurm/evaluate_composition_dependence_v2_cpu.sbatch
EVAL_TEST=tests/unit/test_evaluate_composition_dependence_v2.py
SEALED=${ROOT}/scripts/evaluate_gse267145_axis_count.py
EXT=${ROOT}/scripts/multicovariate_partial.py
PRESPEC_ROOT=${ROOT}/executions/model-check-1060-21183334-composition-prespec-v2
PRESPEC=${PRESPEC_ROOT}/prespecification
STAGE1=${ROOT}/executions/model-data-953-21182514-two-axis-substrate/substrate
COMP=${ROOT}/executions/model-data-954-21183077-composition-substrate/substrate
ARM_B=${ROOT}/executions/model-data-880-21130257-gse135251-source/source
EXTERNAL=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk

TARGET=${ROOT}/executions/model-check-1061-${SLURM_JOB_ID}-composition-dependence-v2
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
bash -n "${WRAPPER}"

# Every file the manifest will name must exist before any work begins.
for required in "${SCRIPT}" "${WRAPPER}" "${EVAL_TEST}" "${SEALED}" "${EXT}"; do
  [[ -f "${required}" ]] || { printf 'Missing: %s\n' "${required}" >&2; exit 1; }
done

cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
344b20137e2436d0736a6bce9f15e6d13849c42d01e43c31f51001f5ccabe0a6  executions/model-check-1060-21183334-composition-prespec-v2/prespecification/ARTIFACTS.json
97f6c12407ac5d093afc978a777983b14804d971963c3291d6510059f9112623  executions/model-data-954-21183077-composition-substrate/substrate/ARTIFACTS.json
68837a8293d5386e6bd7b539dd9107f61d1cfddb783feec373d9853e57be84ab  executions/model-data-953-21182514-two-axis-substrate/substrate/ARTIFACTS.json
5f4e483dc22f1143f1c4447f60e308a5822dd86070540527abbd429e7278e559  scripts/evaluate_gse267145_axis_count.py
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"

# The two held-back v2 trees must still hold byte-identical specification bodies,
# or "duplicate hold back" is the wrong description and this job must not proceed.
# 1057 and 1058 hold one byte-identical body; 1060 supersedes both by carrying
# the two records the evaluator needs. Assert that relationship rather than
# assuming it: the two earlier bodies must still match each other, and the bound
# one must differ from them for the stated reason.
BODY_1057=${ROOT}/executions/model-check-1057-21183286-composition-prespec-v2-v2/prespecification/composition_prespecification_v2.json
BODY_1058=${ROOT}/executions/model-check-1058-21183309-composition-prespec-v2/prespecification/composition_prespecification_v2.json
BODY_BOUND=${PRESPEC}/composition_prespecification_v2.json
cmp -s "${BODY_1057}" "${BODY_1058}" || {
  printf 'The two earlier v2 bodies are not identical; the duplicate-seal record is wrong.\n' >&2
  exit 1
}
if cmp -s "${BODY_1058}" "${BODY_BOUND}"; then
  printf 'The bound body matches the earlier one; the restored keys are missing.\n' >&2
  exit 1
fi
for required_key in closure correction_to_the_frozen_substrate_artifact; do
  grep -q "\"${required_key}\"" "${BODY_BOUND}" || {
    printf 'The bound prespecification lacks %s, which the evaluator reads.\n' "${required_key}" >&2
    exit 1
  }
done

"${PYTHON}" - "${PRESPEC}" "${PRESPEC_ROOT}" "${COMP}" "${STAGE1}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for a in sys.argv[1:]:
    verify_frozen_tree(Path(a).resolve(strict=True))
print("frozen trees verified")
PY

"${PYTHON}" -m unittest -v tests.unit.test_evaluate_composition_dependence_v2 \
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

# The manifest records what actually ran.
sha256sum "${SCRIPT}" "${WRAPPER}" "${EVAL_TEST}" "${SEALED}" "${EXT}" \
  "${PRESPEC}/composition_prespecification_v2.json" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${WRAPPER}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
wrapper = sys.argv[2]
manifest = (root / "source.sha256").read_text(encoding="utf-8")
if wrapper not in manifest:
    raise SystemExit(f"the source manifest does not name the wrapper that ran: {wrapper}")
for stale in ("evaluate_composition_dependence.py\n",
              "evaluate_composition_dependence_cpu.sbatch"):
    if stale in manifest:
        raise SystemExit(f"the source manifest names a v1 sibling that did not run: {stale!r}")

r = root / "composition"
p = json.loads((r / "composition_dependence.json").read_text(encoding="utf-8"))
declared = {"COMPOSITION_INDEPENDENT_COMPONENT_EXISTS",
            "ARMS_DIFFER_WITHIN_RESOLUTION",
            "COMPOSITION_INDEPENDENT_IN_ONE_ARM_ONLY",
            "SET_IS_COMPOSITION_MEDIATED", "INDETERMINATE"}
if p["outcome"] not in declared:
    raise SystemExit(f"undeclared outcome {p['outcome']!r}")
if p["duplicate_seal_record"]["canonical"] != (
        "model-check-1060-21183334-composition-prespec-v2"):
    raise SystemExit("the canonical prespecification is not the one bound")
if not p["the_framing_that_governs_every_reading"][
        "a_composition_dependent_gene_is_not_an_artifact"]:
    raise SystemExit("the mediation framing was dropped")
if len(p["eligible_family"]) != 6:
    raise SystemExit("the family is not the approved six")

# Every arm must contribute, and no arm may vanish silently.
for cohort in ("GSE135251", "GSE130970", "GSE193066"):
    if cohort not in p["arms"]:
        raise SystemExit(f"{cohort} is missing from the result entirely")
    cov = p["arms"][cohort].get("assignment_coverage")
    if cov is None or cov["genes_matching_an_assignment"] == 0:
        raise SystemExit(f"{cohort} matched no assigned gene; a silent arm dropout")
    if cov["background_pool_size"] == 0:
        raise SystemExit(f"{cohort} has an empty background pool")
g = p["criteria"]["t1_a_composition_independent_component_exists"]["gate"]
arms_with_conditions = {c["condition"].split(":")[0] for c in g["conditions"]}
if len(arms_with_conditions) != 3:
    raise SystemExit(f"only {sorted(arms_with_conditions)} produced conditions")
if g["n_conditions_applicable"] == 0 and g["verdict"] != "NO_APPLICABLE_CONDITIONS":
    raise SystemExit("T1 passed with no applicable condition")
if g["passed"] != (g["n_conditions_applicable"] > 0
                   and g["n_conditions_met"] == g["n_conditions_applicable"]):
    raise SystemExit("the T1 verdict does not follow from its conditions")

# The symmetric guard must actually bind: no met condition inside its own MDE.
for c in g["conditions"]:
    if c.get("met") and c.get("margin") is not None and c.get("set_level_mde"):
        if abs(c["margin"]) < c["set_level_mde"]:
            raise SystemExit(
                f"{c['condition']} is met with a margin inside its MDE; the "
                "symmetric guard did not bind")
for cohort, arm in p["arms"].items():
    for axis, cond in arm["conditioning"].items():
        if cond["numerical_rank"] != cond["n_covariates"]:
            raise SystemExit(f"{cohort}/{axis} covariate block is rank deficient")
if p["arms"]["GSE193066"]["n"] != 106:
    raise SystemExit("GSE193066 was not restricted to 106 first-biopsy participants")

verify_frozen_tree(r)
freeze_tree(root, {"artifact_class": "composition_dependence_v2_campaign",
                   "stage": "stage_3_v2", "job": "2_of_2",
                   "outcome": p["outcome"],
                   "prespecification_verified_before_computation": True,
                   "manifest_names_the_wrapper_that_ran": True,
                   "status": "passed"})
verify_frozen_tree(root)
print(json.dumps({"outcome": p["outcome"]}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/composition/ARTIFACTS.json"
