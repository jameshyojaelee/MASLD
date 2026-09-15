#!/usr/bin/env bash
#SBATCH --job-name=model-check-091
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Prove the GSE83452 inductive separation on the GPL16686 code path, before any
# GPL16686 matrix exists.
#
# The GPL570 pass is not evidence here. oligo and SCAN.UPC are different
# packages with different summarization code, and oligo::rma summarizes by
# median polish over the probe-by-array matrix, which borrows across arrays by
# construction. The hazard is assumed present until measured absent, and the
# receipt records gpl570_evidence_carried_over: false so that assumption cannot
# be quietly relaxed.
#
# Pipeline contexts read one CEL per call and must be bit-identical. Multi-array
# estimator contexts are negative controls that must keep diverging. The 1e-12
# perturbation control and the distinct-target-digest check are identical to
# GPL570's.
#
# Raw CEL files are staged onto node-local scratch, checked against the frozen
# raw-source SHA-256 values, and removed on exit.

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
CEL_TAR=${RAW_AUDIT}/raw/GSE83452_RAW.tar
CEL_MANIFEST=${RAW_AUDIT}/GSE83452_cel_manifest.tsv
CEL_MANIFEST_SHA256=a7c20245ed131b27c56b256c450e9337d5f32d51a54e1b12e14f5e78339988f1

PLATFORM_AUDIT=${ROOT}/executions/model-data-067-21080193-platform
PLATFORM_AUDIT_SHA256=0ce6af4196d2fb0905d6ca3d07fbbeb3bf1ce5b37776a63a726cc1d21794f4bf

PREPROCESSING=${ROOT}/config/evaluation/microarray_transfer_preprocessing.toml
PREPROCESSING_SHA256=11162c381017bd16736433f478c0c82d97726ae2d896957d27efa5aee899c27e
LABEL_SEMANTICS=${ROOT}/config/evaluation/microarray_transfer_label_semantics.toml
LABEL_SEMANTICS_SHA256=b66c20ef5c6bc2ce7597affe7816358bd7fd4ca65df8e3f4be1096dbddaf0513

TARGET=${ROOT}/executions/model-check-091-${SLURM_JOB_ID}-firewall
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
CEL_WORK=${TMPDIR:-/tmp}/masld_gse83452_${SLURM_JOB_ID}

cleanup() {
  rm -rf "${CEL_WORK}"
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed GSE83452 firewall stage: %s\n' "${FAILED}" >&2
  fi
}
trap cleanup EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite GSE83452 firewall artifact: %s\n' "${path}" >&2
    exit 1
  fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/validation" "${STAGE}/firewall"
install -d -m 0700 "${CEL_WORK}"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}"         "${RUNTIME}/ARTIFACTS.json" \
  "${RAW_AUDIT_SHA256}"       "${RAW_AUDIT}/ARTIFACTS.json" \
  "${CEL_MANIFEST_SHA256}"    "${CEL_MANIFEST}" \
  "${PLATFORM_AUDIT_SHA256}"  "${PLATFORM_AUDIT}/ARTIFACTS.json" \
  "${PREPROCESSING_SHA256}"   "${PREPROCESSING}" \
  "${LABEL_SEMANTICS_SHA256}" "${LABEL_SEMANTICS}" \
  | sha256sum --check --strict

cd "${ROOT}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

module purge
module load "${PYTHON_MODULE}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 -m unittest -v \
  tests.unit.test_gse83452_oligo_firewall_receipt \
  tests.unit.test_stage_microarray_cel_arrays \
  tests.contract.test_microarray_transfer_contract \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

python3 scripts/stage_microarray_cel_arrays.py \
  --tar "${CEL_TAR}" \
  --manifest "${CEL_MANIFEST}" \
  --output "${CEL_WORK}" \
  --series GSE83452 \
  --expected-records 231 \
  --platform-id GPL16686 \
  --receipt "${STAGE}/firewall/cel_staging_receipt.json" \
  > "${STAGE}/validation/cel_staging_stdout.json" \
  2> "${STAGE}/validation/cel_staging_stderr.txt"

module purge
module load "${R_MODULE}"
R --version > "${STAGE}/validation/R.version.txt"
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/verify_gse83452_oligo_single_array_firewall.R" \
  "${CEL_WORK}" \
  "${CEL_MANIFEST}" \
  "${STAGE}/firewall" \
  > "${STAGE}/validation/oligo_firewall_stdout.txt" \
  2> "${STAGE}/validation/oligo_firewall_stderr.txt"
cp "${STAGE}/firewall/R_sessionInfo.txt" "${STAGE}/validation/R_sessionInfo.txt"

module purge
module load "${PYTHON_MODULE}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 scripts/verify_gse83452_oligo_firewall_receipt.py \
  --receipt "${STAGE}/firewall/oligo_single_array_firewall.json" \
  --output "${STAGE}/firewall/firewall_verdict.json" \
  > "${STAGE}/validation/firewall_verdict_stdout.json" \
  2> "${STAGE}/validation/firewall_verdict_stderr.txt"
python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"

rm -rf "${CEL_WORK}"

sha256sum \
  scripts/stage_microarray_cel_arrays.py \
  scripts/verify_gse83452_oligo_single_array_firewall.R \
  scripts/verify_gse83452_oligo_firewall_receipt.py \
  tests/unit/test_gse83452_oligo_firewall_receipt.py \
  tests/unit/test_stage_microarray_cel_arrays.py \
  slurm/verify_gse83452_oligo_single_array_firewall_cpu.sbatch \
  config/evaluation/microarray_transfer_preprocessing.toml \
  config/evaluation/microarray_transfer_label_semantics.toml \
  "${RUNTIME}/ARTIFACTS.json" \
  "${RAW_AUDIT}/ARTIFACTS.json" \
  "${CEL_MANIFEST}" \
  > "${STAGE}/source.sha256"

if find "${STAGE}" -type l -print -quit | grep -q .; then
  printf 'GSE83452 firewall stage contains a symlink\n' >&2
  exit 1
fi
if find "${STAGE}" -type f -name '*.CEL' -print -quit | grep -q .; then
  printf 'GSE83452 firewall stage contains a raw CEL file\n' >&2
  exit 1
fi

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
    or verdict["series"] != "GSE83452"
    or verdict["platform_id"] != "GPL16686"
    or verdict["method"] != "single_array_core_summary_normalize_false"
    or verdict["held_array_application"] != "one_array_at_a_time_with_frozen_training_object"
    or verdict["targets_verified"] < 3
    or verdict["pipeline_contexts_per_target"] != 4
    or verdict["arrays_per_summarization_call"] != 1
    or verdict["held_array_summary_invariant_to_other_arrays_in_the_run"] is not True
    or verdict["held_array_summary_invariant_to_processing_order"] is not True
    or verdict["perturbation_control_detected_for_every_target"] is not True
    or verdict["distinct_targets_have_distinct_digests"] is not True
    or verdict["multi_array_call_used_as_preprocessing_path"] is not False
    or verdict["gpl570_evidence_carried_over"] is not False
    or verdict["all_sample_RMA_run"] is not False
    or verdict["across_array_quantile_normalization_run"] is not False
    or verdict["GEO_series_matrix_read"] is not False
    or verdict["labels_read"] is not False
    or verdict["model_training_activated"] is not False
    or verdict["sealed_outcomes_read"] is not False
):
    raise SystemExit("GSE83452 single-array firewall verdict differs")
if any(
    record["pipeline_maximum_absolute_difference"] != 0.0
    for record in verdict["per_target"]
):
    raise SystemExit("a target array differs across pipeline contexts")
if staging["staged_arrays"] != 231 or staging["every_array_matches_frozen_sha256"] is not True:
    raise SystemExit("GSE83452 CEL staging receipt differs")

metadata = {
    "artifact_class": "gse83452_gpl16686_single_array_inductive_firewall",
    "series": "GSE83452",
    "platform_id": "GPL16686",
    "cohort_family_id": "antwerp_inserm_shared",
    "method": "single_array_core_summary_normalize_false",
    "held_array_application": "one_array_at_a_time_with_frozen_training_object",
    "core_features": verdict["core_features"],
    "targets_verified": verdict["targets_verified"],
    "pipeline_contexts_per_target": verdict["pipeline_contexts_per_target"],
    "arrays_per_summarization_call": 1,
    "held_array_summary_invariant_to_other_arrays_in_the_run": True,
    "held_array_summary_invariant_to_processing_order": True,
    "multi_array_call_is_batch_dependent": verdict["multi_array_call_is_batch_dependent"],
    "multi_array_call_used_as_preprocessing_path": False,
    "gpl570_evidence_carried_over": False,
    "perturbation_control_detected_for_every_target": True,
    "cel_arrays_verified_against_frozen_sha256": staging["staged_arrays"],
    "scan_single_array_status": (verdict.get("scan_probe") or {}).get("single_array_status"),
    "all_sample_RMA_run": False,
    "across_array_quantile_normalization_run": False,
    "GEO_series_matrix_read": False,
    "labels_read": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
}
freeze_tree(firewall, metadata)
verify_frozen_tree(firewall)
freeze_tree(stage, metadata | {"artifact_class": "gse83452_gpl16686_inductive_firewall_campaign"})
verify_frozen_tree(stage)
PY

trap - EXIT
rm -rf "${CEL_WORK}"
mv "${STAGE}" "${TARGET}"
printf 'GSE83452 GPL16686 single-array firewall: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/firewall/ARTIFACTS.json"
