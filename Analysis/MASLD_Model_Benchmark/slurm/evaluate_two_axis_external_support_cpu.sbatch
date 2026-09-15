#!/usr/bin/env bash
#SBATCH --job-name=model-check-1053
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 1, job 3 of 3: the external-support check.
#
# The prespecification frozen by model-check-1052 (job 21182536) is verified by
# digest before any external outcome association is computed. Both external
# metadata files carry a ragged N+1 header, so every column is resolved by
# measured position and its contents validated against the ordinal it claims.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
SCRIPT=${ROOT}/scripts/evaluate_two_axis_external_support.py
ARM_A=${ROOT}/scripts/evaluate_gse267145_axis_count.py
PRESPEC_ROOT=${ROOT}/executions/model-check-1052-21182536-two-axis-external-prespec
PRESPEC=${PRESPEC_ROOT}/prespecification
SUBSTRATE=${ROOT}/executions/model-data-953-21182514-two-axis-substrate/substrate
EXTERNAL=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk

TARGET=${ROOT}/executions/model-check-1053-${SLURM_JOB_ID}-two-axis-external
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the external-support result.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
bash -n slurm/evaluate_two_axis_external_support_cpu.sbatch

cat > "${STAGE}/frozen_inputs.sha256" <<'SUMS'
b1a687c886c68a4f16efc60ebd4e3c2ecd6c676f6c00508557798cfbba52a5b0  executions/model-check-1052-21182536-two-axis-external-prespec/prespecification/ARTIFACTS.json
68837a8293d5386e6bd7b539dd9107f61d1cfddb783feec373d9853e57be84ab  executions/model-data-953-21182514-two-axis-substrate/substrate/ARTIFACTS.json
5f26597ab6b6aeeaff867d54112554a33d4d0dc35be5f333e17e9010a7199fed  executions/model-check-1050-21182309-arm-b-axis-count/arm_b_axis_count/ARTIFACTS.json
3c5d3cc4e3f47b331c5f8ece04c213c304d49faa3ce3b974cef878c57e4f6ca4  config/evaluation/two_axis_activity_fibrosis_external_activation.json
SUMS
sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"
cat > "${STAGE}/external_inputs.sha256" <<'SUMS'
0d9aa05ce52b7c16501a947e12e25a9b0513a83a2a3a26751307de979a500334  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE130970/GSE130970_counts.tsv
a1b06183e9c9f4444db264afcf2fd18d0215a3d68a7163fd3bc0d3851acd3667  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE193066/GSE193066_counts.tsv
8e8e0661c03a431f7006cf06c1e9bfdf62e4030392c505a4fd93e5cda1a16cba  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE130970/GSE130970_metadata.tsv
667782f577e8ad7cc3d80d5b776d02a4e32c20e8fe54a964ad3102069e1ed7d9  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk/GSE193066/GSE193066_metadata.tsv
SUMS
sha256sum --check --strict "${STAGE}/external_inputs.sha256"

"${PYTHON}" - "${PRESPEC}" "${PRESPEC_ROOT}" "${SUBSTRATE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for argument in sys.argv[1:]:
    verify_frozen_tree(Path(argument).resolve(strict=True))
print("frozen trees verified")
PY

"${PYTHON}" -m unittest -v tests.unit.test_evaluate_two_axis_external_support \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" "${SCRIPT}" \
  --prespecification "${PRESPEC}" \
  --analysis-script "${ARM_A}" \
  --substrate "${SUBSTRATE}" \
  --external-root "${EXTERNAL}" \
  --output "${STAGE}/external" \
  --null-draws 2000 \
  --mde-draws 400 \
  --seed 20260827 \
  > "${STAGE}/external_stdout.txt" 2> "${STAGE}/external_stderr.txt"

sha256sum --check --strict "${STAGE}/frozen_inputs.sha256"
sha256sum --check --strict "${STAGE}/external_inputs.sha256"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum "${SCRIPT}" "${ARM_A}" \
  slurm/evaluate_two_axis_external_support_cpu.sbatch \
  tests/unit/test_evaluate_two_axis_external_support.py \
  "${PRESPEC}/two_axis_external_prespecification.json" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
result = root / "external"
payload = json.loads((result / "two_axis_external_support.json").read_text(encoding="utf-8"))
declared = {"EXTERNALLY_SUPPORTED", "SUPPORTED_SINGLE_ARM",
            "TRANSFERS_WITHOUT_SPECIFICITY", "NOT_EXTERNALLY_SUPPORTED", "INDETERMINATE"}
if payload["outcome"] not in declared:
    raise SystemExit(f"undeclared outcome {payload['outcome']!r}")
if payload["external_support_is_set_level_not_per_gene"] is not True:
    raise SystemExit("the set-level constraint was dropped")
audit = payload["ragged_header_audit"]["per_cohort"]
for cohort, expected in (("GSE130970", 1), ("GSE193066", 1)):
    if audit[cohort]["measured_offset"] != expected:
        raise SystemExit(f"{cohort} offset is not the measured one")
g = payload["arms"]["GSE193066"]["unit_of_inference"]
if g["n"] != 106 or g["distinct_participants_asserted"] != 106:
    raise SystemExit("GSE193066 was not restricted to 106 first-biopsy participants")
if payload["arms"]["GSE130970"]["unit_of_inference"]["distinct_people_assertable"] is not False:
    raise SystemExit("GSE130970 people were made assertable")
if payload["arms"]["GSE193066"]["axes_this_arm_serves"] != ["activity"]:
    raise SystemExit("GSE193066 served an axis it is not scoped for")
for name, gate in ((k, v["gate"]) for k, v in payload["criteria"].items()):
    if gate["n_conditions_applicable"] == 0 and gate["verdict"] != "NO_APPLICABLE_CONDITIONS":
        raise SystemExit(f"{name} passed with no applicable condition")
    if gate["passed"] != (gate["n_conditions_applicable"] > 0
                          and gate["n_conditions_met"] == gate["n_conditions_applicable"]):
        raise SystemExit(f"{name} verdict does not follow from its conditions")
if payload["multiplicity_determination"]["is_a_gate"] is not False:
    raise SystemExit("the multiplicity determination became a gate")
verify_frozen_tree(result)
freeze_tree(root, {"artifact_class": "two_axis_external_support_campaign",
                   "stage": "stage_1", "job": "3_of_3",
                   "outcome": payload["outcome"],
                   "prespecification_verified_before_computation": True,
                   "external_support_is_set_level": True, "status": "passed"})
verify_frozen_tree(root)
print(json.dumps({"outcome": payload["outcome"]}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/external/ARTIFACTS.json"
