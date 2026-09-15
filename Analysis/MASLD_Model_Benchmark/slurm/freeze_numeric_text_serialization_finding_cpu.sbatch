#!/usr/bin/env bash
#SBATCH --job-name=model-check-093
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Freeze the reusable finding that R's default text serialization is lossy for
# doubles, together with the methodological trap that the first verification
# attempt fell into. Both apply to GSE83452's 231 arrays and to every numeric
# matrix this campaign writes as text.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ROUNDTRIP=${ROOT}/executions/model-check-092-21111006-roundtrip/matrix_text_roundtrip.json
ROUNDTRIP_SHA256=ab9afec360e3c9527efc6d848659304e017f7c1ecf4f84b92b6ea9899e452944
SUMMARIZATION=${ROOT}/executions/model-data-086-21110619-summarization
SUMMARIZATION_SHA256=93afd6772e2a5569fdcb9750c49e1549d51d460460b44eb0df67947150a0077a
VERDICT=${SUMMARIZATION}/summary/summarization_verdict.json

TARGET=${ROOT}/executions/model-check-093-${SLURM_JOB_ID}-serialization
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite: %s\n' "${path}" >&2
    exit 1
  fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/validation"

printf '%s  %s\n' \
  "${ROUNDTRIP_SHA256}"     "${ROUNDTRIP}" \
  "${SUMMARIZATION_SHA256}" "${SUMMARIZATION}/summary/ARTIFACTS.json" \
  | sha256sum --check --strict

cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"

python3 -m unittest -v \
  tests.unit.test_build_numeric_text_serialization_finding \
  tests.unit.test_verify_gse49541_summarization_receipt \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

python3 scripts/build_numeric_text_serialization_finding.py \
  --roundtrip-receipt "${ROUNDTRIP}" \
  --summarization-verdict "${VERDICT}" \
  --output "${STAGE}/finding" \
  > "${STAGE}/validation/finding_stdout.json" \
  2> "${STAGE}/validation/finding_stderr.txt"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/build_numeric_text_serialization_finding.py \
  scripts/verify_matrix_text_roundtrip.R \
  scripts/summarize_gse49541_gpl570_single_array.R \
  tests/unit/test_build_numeric_text_serialization_finding.py \
  slurm/freeze_numeric_text_serialization_finding_cpu.sbatch \
  "${ROUNDTRIP}" \
  "${SUMMARIZATION}/summary/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1]).resolve(strict=True)
finding_dir = stage / "finding"
finding = json.loads(
    (finding_dir / "numeric_text_serialization_finding.json").read_text(encoding="utf-8")
)
if (
    finding["status"] != "resolved_by_writing_17_significant_digits_and_asserting_the_round_trip"
    or finding["default_writer_bitwise_identical"] is not False
    or finding["fifteen_digit_writer_bitwise_identical"] is not False
    or finding["seventeen_digit_writer_bitwise_identical"] is not True
    or finding["resolution"]["guard_retained"] is not True
    or finding["fraction_differing_under_default_writer"] < 0.5
    or "methodological_trap" not in finding
):
    raise SystemExit("serialization finding differs from the bound evidence")

metadata = {
    "artifact_class": "numeric_text_serialization_finding",
    "finding_id": finding["finding_id"],
    "values_tested": finding["values_tested"],
    "fraction_differing_under_default_writer": finding[
        "fraction_differing_under_default_writer"
    ],
    "maximum_absolute_error_log2": finding["maximum_absolute_error_log2"],
    "resolution": "sprintf_17g_plus_bitwise_round_trip_assertion",
    "guard_retained": True,
    "applies_to": finding["generalization"]["named_exposures"],
    "labels_read": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
}
freeze_tree(finding_dir, metadata)
verify_frozen_tree(finding_dir)
freeze_tree(stage, metadata | {"artifact_class": "numeric_text_serialization_finding_campaign"})
verify_frozen_tree(stage)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Numeric text serialization finding: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/finding/ARTIFACTS.json"
