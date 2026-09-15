#!/usr/bin/env bash
#SBATCH --job-name=model-check-089
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Prove the GSE49541 inductive separation on real CEL files before any expression
# matrix is produced. Re-run after model-check-085 (job 21107058) failed.
#
# That failure stands and is not weakened here. It showed that a multi-array
# frma() call is batch dependent, because the default robust_weighted_average
# re-estimates its robustness weights over every array in the call. The fixture
# now separates the two questions that were previously conflated.
#
# Pipeline contexts read exactly ONE CEL per frma() call, which is what the
# frozen requirements require, and differ only in what else the run processed and
# in what order. Those must be bit-identical: that is the separation.
#
# Estimator contexts keep the original multi-array calls as negative controls.
# They are expected to diverge, and the verdict module REQUIRES them to keep
# diverging, so batching cannot be reintroduced without the fixture noticing.
#
# The 1e-12 perturbation control is carried over unchanged, and distinct targets
# must still carry distinct digests.
#
# Raw CEL files are staged onto node-local scratch, checked against the frozen
# raw-source SHA-256 values, and removed on exit. They never enter a published
# output-file tree, so the source's redistribution prohibition holds.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export R_LIBS_USER=/nonexistent/masld_bench_r_user_library

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT_ROOT}/Analysis/MASLD_Model_Benchmark

RUNTIME=${ROOT}/executions/environments/microarray-r-bioc-21104442
RUNTIME_SHA256=8115e0966804dfc4b7bfac3f95016896d099d82fc51f3428fab1ec34182a7d07
R_MODULE=R/4.5.2-mba
PYTHON_MODULE=python/3.11.5-GCCcore-13.2.0

RAW_AUDIT=${ROOT}/executions/model-data-065-21079931
RAW_AUDIT_SHA256=8b30a2265c421f2917ed886d99cfaa08ee8d874abd40e68c8c72d48c726781d1
CEL_TAR=${RAW_AUDIT}/raw/GSE49541_RAW.tar
CEL_MANIFEST=${RAW_AUDIT}/GSE49541_cel_manifest.tsv
CEL_MANIFEST_SHA256=46204986a7bbefbcde169411d7fb5b6bded4b6babd560e49fc3e47e065c579f8

PLATFORM_AUDIT=${ROOT}/executions/model-data-067-21080193-platform
PLATFORM_AUDIT_SHA256=0ce6af4196d2fb0905d6ca3d07fbbeb3bf1ce5b37776a63a726cc1d21794f4bf
PLATFORM_AXIS=${PLATFORM_AUDIT}/GPL570_selected_annotation.tsv
PLATFORM_AXIS_SHA256=4bbd91765d1a80b12439e8c39b2aca9218a4135a80fc6a50f6f4eed276d61b69

PREPROCESSING=${ROOT}/config/evaluation/microarray_transfer_preprocessing.toml
PREPROCESSING_SHA256=11162c381017bd16736433f478c0c82d97726ae2d896957d27efa5aee899c27e
TASKSPEC=${ROOT}/config/evaluation/gse49541_fibrosis_transfer_task.toml
TASKSPEC_SHA256=9b11584b3c1b2a01d4bb262790274032c31f92ff0b1da02a66588aa3e75a6ec3

TARGET=${ROOT}/executions/model-check-089-${SLURM_JOB_ID}-firewall
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
CEL_WORK=${TMPDIR:-/tmp}/masld_gse49541_${SLURM_JOB_ID}

cleanup() {
  rm -rf "${CEL_WORK}"
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed GSE49541 firewall stage: %s\n' "${FAILED}" >&2
  fi
}
trap cleanup EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite GSE49541 firewall artifact: %s\n' "${path}" >&2
    exit 1
  fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/validation" "${STAGE}/firewall"
install -d -m 0700 "${CEL_WORK}"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}"        "${RUNTIME}/ARTIFACTS.json" \
  "${RAW_AUDIT_SHA256}"      "${RAW_AUDIT}/ARTIFACTS.json" \
  "${CEL_MANIFEST_SHA256}"   "${CEL_MANIFEST}" \
  "${PLATFORM_AUDIT_SHA256}" "${PLATFORM_AUDIT}/ARTIFACTS.json" \
  "${PLATFORM_AXIS_SHA256}"  "${PLATFORM_AXIS}" \
  "${PREPROCESSING_SHA256}"  "${PREPROCESSING}" \
  "${TASKSPEC_SHA256}"       "${TASKSPEC}" \
  | sha256sum --check --strict

cd "${ROOT}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

## ------------------------------------------------------------------ unit tests
module purge
module load "${PYTHON_MODULE}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 -m unittest -v \
  tests.unit.test_gse49541_frma_firewall_receipt \
  tests.unit.test_stage_gse49541_cel_arrays \
  tests.contract.test_microarray_transfer_contract \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

## -------------------------------- stage the whole cohort onto node-local disk
python3 scripts/stage_gse49541_cel_arrays.py \
  --tar "${CEL_TAR}" \
  --manifest "${CEL_MANIFEST}" \
  --output "${CEL_WORK}" \
  --receipt "${STAGE}/firewall/cel_staging_receipt.json" \
  > "${STAGE}/validation/cel_staging_stdout.json" \
  2> "${STAGE}/validation/cel_staging_stderr.txt"

## ------------------------------ single-array isolation on the real GPL570 CELs
module purge
module load "${R_MODULE}"
R --version > "${STAGE}/validation/R.version.txt"
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/verify_gse49541_frma_single_array_firewall.R" \
  "${CEL_WORK}" \
  "${CEL_MANIFEST}" \
  "${PLATFORM_AXIS}" \
  "${STAGE}/firewall" \
  > "${STAGE}/validation/frma_firewall_stdout.txt" \
  2> "${STAGE}/validation/frma_firewall_stderr.txt"
cp "${STAGE}/firewall/R_sessionInfo.txt" "${STAGE}/validation/R_sessionInfo.txt"

## ------------------------------------------------------------------- verdict
module purge
module load "${PYTHON_MODULE}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 scripts/verify_gse49541_frma_firewall_receipt.py \
  --receipt "${STAGE}/firewall/frma_single_array_firewall.json" \
  --output "${STAGE}/firewall/firewall_verdict.json" \
  > "${STAGE}/validation/firewall_verdict_stdout.json" \
  2> "${STAGE}/validation/firewall_verdict_stderr.txt"
python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"

rm -rf "${CEL_WORK}"

sha256sum \
  scripts/stage_gse49541_cel_arrays.py \
  scripts/verify_gse49541_frma_single_array_firewall.R \
  scripts/verify_gse49541_frma_firewall_receipt.py \
  tests/unit/test_gse49541_frma_firewall_receipt.py \
  tests/unit/test_stage_gse49541_cel_arrays.py \
  slurm/verify_gse49541_frma_single_array_firewall_cpu.sbatch \
  config/evaluation/microarray_transfer_preprocessing.toml \
  config/evaluation/gse49541_fibrosis_transfer_task.toml \
  "${RUNTIME}/ARTIFACTS.json" \
  "${RAW_AUDIT}/ARTIFACTS.json" \
  "${CEL_MANIFEST}" \
  "${PLATFORM_AXIS}" \
  > "${STAGE}/source.sha256"

if find "${STAGE}" -type l -print -quit | grep -q .; then
  printf 'GSE49541 firewall stage contains a symlink\n' >&2
  find "${STAGE}" -type l -print >&2
  exit 1
fi
if find "${STAGE}" -type f -name '*.CEL' -print -quit | grep -q .; then
  printf 'GSE49541 firewall stage contains a raw CEL file\n' >&2
  exit 1
fi

## ---------------------------------------------------------------------- freeze
python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1]).resolve(strict=True)
firewall = stage / "firewall"
verdict = json.loads((firewall / "firewall_verdict.json").read_text(encoding="utf-8"))
staging = json.loads((firewall / "cel_staging_receipt.json").read_text(encoding="utf-8"))
if (
    verdict["status"] != "pass_inductive_firewall_on_real_CEL"
    or verdict["series"] != "GSE49541"
    or verdict["platform_id"] != "GPL570"
    or verdict["method"] != "frma_with_exact_hgu133plus2frmavecs"
    or verdict["held_array_application"] != "one_array_at_a_time_with_frozen_training_object"
    or verdict["targets_verified"] < 3
    or verdict["pipeline_contexts_per_target"] != 4
    or verdict["arrays_per_frma_call"] != 1
    or verdict["held_array_summary_invariant_to_other_arrays_in_the_run"] is not True
    or verdict["held_array_summary_invariant_to_processing_order"] is not True
    or verdict["perturbation_control_detected_for_every_target"] is not True
    or verdict["distinct_targets_have_distinct_digests"] is not True
    or verdict["multi_array_frma_call_is_batch_dependent"] is not True
    or verdict["multi_array_frma_call_used_as_preprocessing_path"] is not False
    or verdict["all_sample_RMA_run"] is not False
    or verdict["across_array_quantile_normalization_run"] is not False
    or verdict["GEO_series_matrix_read"] is not False
    or verdict["labels_read"] is not False
    or verdict["model_training_activated"] is not False
    or verdict["sealed_outcomes_read"] is not False
):
    raise SystemExit("GSE49541 single-array firewall verdict differs")
if any(
    record["pipeline_maximum_absolute_difference"] != 0.0
    for record in verdict["per_target"]
):
    raise SystemExit("a target array differs across pipeline contexts")
if all(
    record["estimator_maximum_differing_elements"] == 0
    for record in verdict["per_target"]
):
    raise SystemExit("the multi-array negative control stopped diverging")
if (
    staging["staged_arrays"] != 72
    or staging["every_array_matches_frozen_sha256"] is not True
    or staging["raw_CEL_written_into_published_artifact"] is not False
):
    raise SystemExit("GSE49541 CEL staging receipt differs")
freeze_tree(
    firewall,
    {
        "artifact_class": "gse49541_gpl570_single_array_inductive_firewall",
        "series": "GSE49541",
        "platform_id": "GPL570",
        "cohort_family_id": "gse31803_gse49541_fibrosis_array",
        "method": "frma_with_exact_hgu133plus2frmavecs",
        "held_array_application": "one_array_at_a_time_with_frozen_training_object",
        "targets_verified": verdict["targets_verified"],
        "pipeline_contexts_per_target": verdict["pipeline_contexts_per_target"],
        "arrays_per_frma_call": 1,
        "held_array_summary_invariant_to_other_arrays_in_the_run": True,
        "held_array_summary_invariant_to_processing_order": True,
        "multi_array_frma_call_is_batch_dependent": True,
        "multi_array_frma_call_used_as_preprocessing_path": False,
        "batching_would_break_the_firewall": True,
        "perturbation_control_detected_for_every_target": True,
        "cel_arrays_verified_against_frozen_sha256": staging["staged_arrays"],
        "all_sample_RMA_run": False,
        "across_array_quantile_normalization_run": False,
        "GEO_series_matrix_read": False,
        "labels_read": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(firewall)

freeze_tree(
    stage,
    {
        "artifact_class": "gse49541_gpl570_inductive_firewall_campaign",
        "series": "GSE49541",
        "platform_id": "GPL570",
        "single_array_isolation_proven_on_real_CEL": True,
        "multi_array_frma_call_is_batch_dependent": True,
        "targets_verified": verdict["targets_verified"],
        "labels_read": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(stage)
PY

trap - EXIT
rm -rf "${CEL_WORK}"
mv "${STAGE}" "${TARGET}"
printf 'GSE49541 GPL570 single-array firewall: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/firewall/ARTIFACTS.json"
