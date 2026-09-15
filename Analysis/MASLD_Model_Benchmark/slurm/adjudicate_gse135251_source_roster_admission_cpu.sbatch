#!/usr/bin/env bash
#SBATCH --job-name=model-check-892
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
FROZEN_SOURCE=${ROOT}/executions/model-data-880-21130257-gse135251-source/source
FROZEN_SOURCE_ARTIFACTS_SHA256=ad31235cc6eabc0db2ce3a760d8e0340e9f1d9690f672e79271ac4c3406fa63b
LABEL_CONTRACT=${ROOT}/config/evaluation/microarray_transfer_label_semantics.toml
LABEL_CONTRACT_SHA256=b66c20ef5c6bc2ce7597affe7816358bd7fd4ca65df8e3f4be1096dbddaf0513
TASKSPEC=${ROOT}/config/evaluation/gse83452_baseline_nash_transfer_task.toml
TASKSPEC_SHA256=c567022d32d478b364d6a3cd9c84ba692a023bac930d99d3bb262499d60e943b
USER_APPROVAL_DATE=2026-08-26
TARGET=${ROOT}/executions/model-check-892-${SLURM_JOB_ID}-gse135251-admission
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the GSE135251 admission record.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FROZEN_SOURCE_ARTIFACTS_SHA256}" "${FROZEN_SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/adjudicate_gse135251_source_roster_admission_cpu.sbatch

# The two frozen bases, BEFORE.  Captured to a file so the after-check is a
# comparison against evidence on disk, not against a value re-typed here.
sha256sum "${LABEL_CONTRACT}" "${TASKSPEC}" > "${STAGE}/frozen_bases_before.sha256"
printf '%s  %s\n' "${LABEL_CONTRACT_SHA256}" "${LABEL_CONTRACT}" | sha256sum --check --strict
printf '%s  %s\n' "${TASKSPEC_SHA256}" "${TASKSPEC}" | sha256sum --check --strict

# The adjudicator may never write to either frozen base.
if grep -nE '"[wax]"|open\(.*mode=.*w|write_text|\.unlink|shutil\.(copy|move)' \
    scripts/adjudicate_gse135251_source_roster_admission.py \
    | grep -vE 'output|handle\.write|write_tsv|"x"'; then
  printf 'The adjudicator contains an unexpected write path.\n' >&2
  exit 1
fi

"${PYTHON}" - "${FROZEN_SOURCE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
for name in ("", "raw", "molecular", "outcomes", "folds"):
    verify_frozen_tree(root / name if name else root)
print("frozen GSE135251 source verification passed")
PY

"${PYTHON}" -m unittest -v \
  tests.unit.test_adjudicate_gse135251_source_roster_admission \
  > "${STAGE}/unit_tests.txt" 2>&1

"${PYTHON}" scripts/adjudicate_gse135251_source_roster_admission.py \
  --frozen-source "${FROZEN_SOURCE}" \
  --label-contract "${LABEL_CONTRACT}" \
  --taskspec "${TASKSPEC}" \
  --user-approval-date "${USER_APPROVAL_DATE}" \
  --output "${STAGE}/admission" \
  > "${STAGE}/adjudication_stdout.json" 2> "${STAGE}/adjudication_stderr.txt"

# The two frozen bases, AFTER.  Byte-identical or this job fails.
sha256sum "${LABEL_CONTRACT}" "${TASKSPEC}" > "${STAGE}/frozen_bases_after.sha256"
if ! cmp -s "${STAGE}/frozen_bases_before.sha256" "${STAGE}/frozen_bases_after.sha256"; then
  printf 'A frozen base changed during adjudication.\n' >&2
  diff "${STAGE}/frozen_bases_before.sha256" "${STAGE}/frozen_bases_after.sha256" >&2 || true
  exit 1
fi
sha256sum --check --strict "${STAGE}/frozen_bases_after.sha256"

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/adjudicate_gse135251_source_roster_admission.py \
  tests/unit/test_adjudicate_gse135251_source_roster_admission.py \
  slurm/adjudicate_gse135251_source_roster_admission_cpu.sbatch \
  "${LABEL_CONTRACT}" \
  "${TASKSPEC}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${FROZEN_SOURCE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
admission = root / "admission"
record = json.loads((admission / "admission_record.json").read_text(encoding="utf-8"))
with (admission / "nas_threshold_separability.tsv").open(encoding="utf-8", newline="") as handle:
    grid = list(csv.DictReader(handle, delimiter="\t"))

item1 = record["evidence"]["item_1_nafl_call_is_not_a_nas_threshold"]
clause = record["clause_adjudication"]
if (
    record.get("status") != "admitted_to_the_gse83452_source_roster"
    or record.get("record_is_additive_overlay") is not True
    or record.get("frozen_bases_edited") is not False
    or record["admitted_contrast"]["positive_n"] != 131
    or record["admitted_contrast"]["negative_n"] != 41
    or record["admitted_contrast"]["excluded_source_groups"] != ["control"]
    or item1["call_is_a_nas_sum_threshold"] is not False
    or item1["best_accuracy"] >= 1.0
    or clause["negative_required"]["met"] is not True
    or clause["negative_required"]["basis"] != "recomputed_on_disk_evidence"
    or clause["positive_required"]["met_on_disk"] is not False
    or clause["positive_required"]["basis"] != "authoritatively_documented_equivalent_definition"
    or record["evidence"]["item_3_registry_characterization"]["primary_methods_text_read"] is not False
    or record["age_sex_logistic"]["fittable_from_this_source"] is not False
    or record["admission"]["prior_artifact_stamps_rewritten"] is not False
    or record["claim_boundary"]["champion_eligible"] is not False
    or record["claim_boundary"]["external_development_only"] is not True
    or len(record["not_established"]) < 4
    or len(grid) != 10
):
    raise SystemExit("GSE135251 admission-record verification differs")

freeze_tree(
    root,
    {
        "artifact_class": "gse135251_source_roster_admission_campaign",
        "source_roster_admission_state": "admitted",
        "frozen_bases_edited": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
verify_frozen_tree(admission)
print(json.dumps({
    "status": record["status"],
    "admitted_contrast": record["admitted_contrast"]["arm_id"],
    "positive_n": record["admitted_contrast"]["positive_n"],
    "negative_n": record["admitted_contrast"]["negative_n"],
    "best_nas_threshold_accuracy": item1["best_accuracy"],
    "best_nas_threshold_misclassified": item1["best_misclassified"],
    "overlapping_nas_values": item1["overlapping_nas_values"],
    "negative_required_basis": clause["negative_required"]["basis"],
    "positive_required_basis": clause["positive_required"]["basis"],
    "positive_required_met_on_disk": clause["positive_required"]["met_on_disk"],
}, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/admission/ARTIFACTS.json"
cat "${TARGET}/frozen_bases_after.sha256"
