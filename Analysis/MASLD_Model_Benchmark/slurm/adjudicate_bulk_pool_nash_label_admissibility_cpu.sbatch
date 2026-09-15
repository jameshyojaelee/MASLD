#!/usr/bin/env bash
#SBATCH --job-name=model-check-893
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
METADATA_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/bulk
LABEL_CONTRACT=${ROOT}/config/evaluation/microarray_transfer_label_semantics.toml
LABEL_CONTRACT_SHA256=b66c20ef5c6bc2ce7597affe7816358bd7fd4ca65df8e3f4be1096dbddaf0513
TASKSPEC=${ROOT}/config/evaluation/gse83452_baseline_nash_transfer_task.toml
TASKSPEC_SHA256=c567022d32d478b364d6a3cd9c84ba692a023bac930d99d3bb262499d60e943b
TARGET=${ROOT}/executions/model-check-893-${SLURM_JOB_ID}-bulk-pool-label-admissibility
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the label-admissibility record.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/adjudicate_bulk_pool_nash_label_admissibility_cpu.sbatch

# The two frozen bases, BEFORE.
sha256sum "${LABEL_CONTRACT}" "${TASKSPEC}" > "${STAGE}/frozen_bases_before.sha256"
printf '%s  %s\n' "${LABEL_CONTRACT_SHA256}" "${LABEL_CONTRACT}" | sha256sum --check --strict
printf '%s  %s\n' "${TASKSPEC_SHA256}" "${TASKSPEC}" | sha256sum --check --strict

# The read-only cohort tree must not be written.  Recorded before and after.
find "${METADATA_ROOT}"/GSE162694 "${METADATA_ROOT}"/GSE240729 "${METADATA_ROOT}"/GSE213621 \
  -maxdepth 1 -type f -name '*_metadata.tsv' -printf '%p %s %T@\n' | sort \
  > "${STAGE}/source_tree_before.txt"

"${PYTHON}" -m unittest -v \
  tests.unit.test_adjudicate_bulk_pool_nash_label_admissibility \
  > "${STAGE}/unit_tests.txt" 2>&1

"${PYTHON}" scripts/adjudicate_bulk_pool_nash_label_admissibility.py \
  --metadata-root "${METADATA_ROOT}" \
  --label-contract "${LABEL_CONTRACT}" \
  --taskspec "${TASKSPEC}" \
  --output "${STAGE}/adjudication" \
  > "${STAGE}/adjudication_stdout.json" 2> "${STAGE}/adjudication_stderr.txt"

# The two frozen bases, AFTER.  Byte-identical or this job fails.
sha256sum "${LABEL_CONTRACT}" "${TASKSPEC}" > "${STAGE}/frozen_bases_after.sha256"
if ! cmp -s "${STAGE}/frozen_bases_before.sha256" "${STAGE}/frozen_bases_after.sha256"; then
  printf 'A frozen base changed during adjudication.\n' >&2
  exit 1
fi
sha256sum --check --strict "${STAGE}/frozen_bases_after.sha256"
find "${METADATA_ROOT}"/GSE162694 "${METADATA_ROOT}"/GSE240729 "${METADATA_ROOT}"/GSE213621 \
  -maxdepth 1 -type f -name '*_metadata.tsv' -printf '%p %s %T@\n' | sort \
  > "${STAGE}/source_tree_after.txt"
if ! cmp -s "${STAGE}/source_tree_before.txt" "${STAGE}/source_tree_after.txt"; then
  printf 'The read-only cohort tree changed during adjudication.\n' >&2
  exit 1
fi

"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/adjudicate_bulk_pool_nash_label_admissibility.py \
  tests/unit/test_adjudicate_bulk_pool_nash_label_admissibility.py \
  slurm/adjudicate_bulk_pool_nash_label_admissibility_cpu.sbatch \
  "${LABEL_CONTRACT}" \
  "${TASKSPEC}" \
  "${RUNTIME}/ARTIFACTS.json" \
  "${METADATA_ROOT}/GSE162694/GSE162694_metadata.tsv" \
  "${METADATA_ROOT}/GSE240729/GSE240729_metadata.tsv" \
  "${METADATA_ROOT}/GSE213621/GSE213621_metadata.tsv" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
adjudication = root / "adjudication"
record = json.loads((adjudication / "admissibility_record.json").read_text(encoding="utf-8"))
with (adjudication / "clause_failures.tsv").open(encoding="utf-8", newline="") as handle:
    failures = list(csv.DictReader(handle, delimiter="\t"))

by_cohort = {row["cohort"]: row for row in record["verdicts"]}
if (
    record.get("status") != "no_cohort_admitted"
    or record.get("frozen_bases_edited") is not False
    or record.get("record_is_additive_overlay") is not True
    or record.get("admitted") != []
    or sorted(record.get("does_not_pass", [])) != ["GSE162694", "GSE213621", "GSE240729"]
    or record["admission"]["source_roster_admission_state"] != "proposed_not_approved"
    or record["admission"]["approved_by"] is not None
    or record["admission"]["prior_artifact_stamps_rewritten"] is not False
    or record["scope_limits"]["a_frozen_fibrosis_label_semantics_contract_exists"] is not False
    or record["scope_limits"]["fibrosis_label_admissibility_adjudicated_here"] is not False
    or by_cohort["GSE162694"]["verdict"] != "does_not_pass"
    or by_cohort["GSE240729"]["verdict"] != "does_not_pass"
    or by_cohort["GSE213621"]["verdict"] != "does_not_pass"
    or by_cohort["GSE162694"]["has_explicit_nash_positive"] is not True
    or by_cohort["GSE240729"]["label_is_fibrosis_stage_only"] is not True
    or by_cohort["GSE213621"]["clauses_failed"] != 2
    or not failures
):
    raise SystemExit("bulk-pool label-admissibility verification differs")

freeze_tree(
    root,
    {
        "artifact_class": "bulk_existing_source_pool_label_admissibility_campaign",
        "cohorts_adjudicated": len(record["cohorts_adjudicated"]),
        "admitted": 0,
        "frozen_bases_edited": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
verify_frozen_tree(adjudication)
print(json.dumps({
    "status": record["status"],
    "admitted": record["admitted"],
    "does_not_pass": record["does_not_pass"],
    "verdicts": {
        c: {
            "verdict": v["verdict"],
            "clauses_failed": v["clauses_failed"],
            "failed": [f["prohibited_substitution"] for f in v["failed_clauses"]],
        }
        for c, v in by_cohort.items()
    },
}, indent=2, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/adjudication/ARTIFACTS.json"
cat "${TARGET}/frozen_bases_after.sha256"
