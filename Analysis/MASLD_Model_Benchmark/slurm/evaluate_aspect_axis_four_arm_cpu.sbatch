#!/usr/bin/env bash
#SBATCH --job-name=aspect-axis-four-arm
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Stage 1 of the four-arm aspect decomposition.
#
# PRESPEC_SHA256 is the digest of the Stage 0 prespecification file. It is
# pinned here by hand AFTER Stage 0 completed and BEFORE any association was
# computed. The evaluation exits without opening a matrix if it does not match.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
PRESPEC_DIR=${ROOT}/executions/aspect-axis-four-arm-prespec-21324227/prespecification
PRESPEC_SHA256=70cb74fd1b06442570d8e841c40a4100bc9a4164b3f531d47e52123581baa128

TARGET=${ROOT}/executions/aspect-axis-four-arm-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the four-arm result.\n' >&2; exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}/scripts:${ROOT}"
bash -n slurm/evaluate_aspect_axis_four_arm_cpu.sbatch

"${PYTHON}" -m unittest -v tests.unit.test_aspect_axis_instrument \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

# G_PRESPEC_DIGEST, demonstrated live: a one-byte tamper on a COPY must change
# the digest and must make the evaluator refuse. If this block prints anything
# other than the refusal, the digest gate is decorative.
TAMPER=${STAGE}/validation/tamper
install -d -m 0750 "${TAMPER}"
cp -r "${PRESPEC_DIR}/." "${TAMPER}/"
printf ' ' >> "${TAMPER}/aspect_axis_prespecification.json"
TAMPERED_SHA=$(sha256sum "${TAMPER}/aspect_axis_prespecification.json" | cut -d' ' -f1)
if [[ "${TAMPERED_SHA}" == "${PRESPEC_SHA256}" ]]; then
  printf 'G_PRESPEC_DIGEST is dead: a tamper did not change the bytes.\n' >&2; exit 1
fi
set +e
"${PYTHON}" scripts/evaluate_aspect_axis_four_arm.py \
  --prespecification "${TAMPER}" \
  --prespecification-sha256 "${PRESPEC_SHA256}" \
  --output "${STAGE}/validation/tamper_output" \
  > "${STAGE}/validation/tamper_refusal.txt" 2>&1
TAMPER_RC=$?
set -e
if [[ ${TAMPER_RC} -eq 0 ]]; then
  printf 'G_PRESPEC_DIGEST is dead: the evaluator accepted a tampered prespecification.\n' >&2
  exit 1
fi
if [[ -e "${STAGE}/validation/tamper_output" ]]; then
  printf 'G_PRESPEC_DIGEST is dead: the evaluator produced output from a tampered prespecification.\n' >&2
  exit 1
fi
printf 'G_PRESPEC_DIGEST fired on a deliberate tamper (rc=%s)\n' "${TAMPER_RC}" \
  | tee -a "${STAGE}/validation/tamper_refusal.txt"
rm -r "${TAMPER}"

find "${PRESPEC_DIR}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_before.txt"

"${PYTHON}" scripts/evaluate_aspect_axis_four_arm.py \
  --prespecification "${PRESPEC_DIR}" \
  --prespecification-sha256 "${PRESPEC_SHA256}" \
  --output "${STAGE}/four_arm" \
  --n-perm 5000 --seed 20260901 \
  > "${STAGE}/four_arm_stdout.txt" 2> "${STAGE}/four_arm_stderr.txt"

find "${PRESPEC_DIR}" -type f -printf '%p %s\n' | sort > "${STAGE}/frozen_after.txt"
cmp -s "${STAGE}/frozen_before.txt" "${STAGE}/frozen_after.txt" \
  || { printf 'The sealed prespecification changed during the run.\n' >&2; exit 1; }

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum scripts/evaluate_aspect_axis_four_arm.py scripts/aspect_axis_common.py \
  scripts/aspect_axis_arms.py tests/unit/test_aspect_axis_instrument.py \
  slurm/evaluate_aspect_axis_four_arm_cpu.sbatch \
  "${PRESPEC_DIR}/aspect_axis_prespecification.json" > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" "${PRESPEC_SHA256}" <<'PY'
import json, pathlib, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = pathlib.Path(sys.argv[1]).resolve(strict=True)
payload = json.loads((root / "four_arm" / "aspect_axis_result.json").read_text())
if payload["prespecification_sha256"] != sys.argv[2]:
    raise SystemExit("the result is not against the pinned prespecification")
declared = {"THREE_AXES", "TWO_AXES", "ONE_AXIS", "NO_AXIS_RESOLVED",
            "INDETERMINATE"}
for arm, s in payload["per_arm"].items():
    if s["three_aspect_outcome"] not in declared:
        raise SystemExit(f"arm {arm}: undeclared outcome {s['three_aspect_outcome']}")
    n1 = s["controls"]["negative_N1_duplicate"]
    if n1 is None or not n1["passes"]:
        raise SystemExit(f"arm {arm}: the duplicate-aspect negative control did "
                         "not return NOT_APPLICABLE_COLLINEAR")
    n2 = s["controls"]["negative_N2_shuffled"]
    if n2 is None or not n2["passes"]:
        raise SystemExit(f"arm {arm}: a shuffled exposure was called UNIQUE")
    for aspect, status in s["aspect_status"].items():
        if status not in {"supported", "indeterminate", "untestable",
                          "not_independent"}:
            raise SystemExit(f"arm {arm}/{aspect}: undeclared status {status}")
verify_frozen_tree(root / "four_arm")
freeze_tree(root, {"artifact_class": "aspect_axis_count_four_arm",
                   "stage": "stage_1", "status": "passed",
                   "prespecification_sha256": sys.argv[2]})
verify_frozen_tree(root)
print(json.dumps({a: {"outcome": s["three_aspect_outcome"],
                      "axes": s["named_aspect_axes_supported"],
                      "status": s["aspect_status"]}
                  for a, s in payload["per_arm"].items()}, indent=2))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
echo "RESULT_DIR=${TARGET}"
