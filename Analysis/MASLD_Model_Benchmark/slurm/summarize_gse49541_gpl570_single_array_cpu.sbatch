#!/usr/bin/env bash
#SBATCH --job-name=model-data-086
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Label-blind single-array fRMA summarization of all 72 GSE49541 GPL570 arrays.
#
# Runs only after the inductive separation has been proven on these same CEL files
# by model-check-085; that output file is bound by SHA-256 below, so this job cannot
# produce a matrix whose single-array property is unproven.
#
# Every array is summarized in its own fRMA call. The gene matrix keeps only the
# frozen one-to-one GENCODE v49 features from model-data-083 and collapses them
# by the training-frozen within-array median. Participant masks are built through
# a column allowlist that cannot see the deposited fibrosis group, so this whole
# job is label-blind by construction rather than by promise.
#
# This job fits no model, opens no label, reads no GEO series matrix, runs no
# all-sample RMA, and runs no across-array normalization.

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

FIREWALL=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-check-089-21107882-firewall/firewall
FIREWALL_SHA256=c292dadc8454f2d13c682205a3dbf2e738082e591fa10ff0d6a9c4806f3a0e29

CROSSWALK=${ROOT}/executions/model-data-083-21106841-crosswalk/crosswalk
CROSSWALK_SHA256=79b93b1f23b895d77cba3df56c68b6afcb85407130383740fa364c0c31483d18
CROSSWALK_TSV=${CROSSWALK}/gpl570_gencode_v49_crosswalk.tsv
CROSSWALK_TSV_SHA256=1aea919a955136f1cd3b01e244e57878e1e18db745449c33295af9c085bc35a1

RAW_AUDIT=${ROOT}/executions/model-data-065-21079931
RAW_AUDIT_SHA256=8b30a2265c421f2917ed886d99cfaa08ee8d874abd40e68c8c72d48c726781d1
CEL_TAR=${RAW_AUDIT}/raw/GSE49541_RAW.tar
CEL_MANIFEST=${RAW_AUDIT}/GSE49541_cel_manifest.tsv
CEL_MANIFEST_SHA256=46204986a7bbefbcde169411d7fb5b6bded4b6babd560e49fc3e47e065c579f8

PLATFORM_AUDIT=${ROOT}/executions/model-data-067-21080193-platform
PLATFORM_AUDIT_SHA256=0ce6af4196d2fb0905d6ca3d07fbbeb3bf1ce5b37776a63a726cc1d21794f4bf
PLATFORM_AXIS=${PLATFORM_AUDIT}/GPL570_selected_annotation.tsv
PLATFORM_AXIS_SHA256=4bbd91765d1a80b12439e8c39b2aca9218a4135a80fc6a50f6f4eed276d61b69

PARTICIPANT_AUDIT=${ROOT}/executions/model-data-067-21080193-participants
PARTICIPANT_AUDIT_SHA256=1ba06484d9d064198c6c85d480c6ae2c0dadea712d0b53f9f15d81c6aa35fddb
PARTICIPANTS=${PARTICIPANT_AUDIT}/gse49541_participants.tsv
PARTICIPANTS_SHA256=da8985e1988e49fe4804c6c961719feeaf8cab06dd466a7fbc340accf9f9df27

PREPROCESSING=${ROOT}/config/evaluation/microarray_transfer_preprocessing.toml
PREPROCESSING_SHA256=11162c381017bd16736433f478c0c82d97726ae2d896957d27efa5aee899c27e
TASKSPEC=${ROOT}/config/evaluation/gse49541_fibrosis_transfer_task.toml
TASKSPEC_SHA256=9b11584b3c1b2a01d4bb262790274032c31f92ff0b1da02a66588aa3e75a6ec3
GATE=${ROOT}/config/evaluation/gse49541_fibrosis_transfer_promotion_gate.json
GATE_SHA256=791dfa3199c3bdbdeb74bfae429cc8820364862c7187f5d1dbe11b46c03dc384

TARGET=${ROOT}/executions/model-data-086-${SLURM_JOB_ID}-summarization
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
CEL_WORK=${TMPDIR:-/tmp}/masld_gse49541_${SLURM_JOB_ID}

cleanup() {
  rm -rf "${CEL_WORK}"
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed GSE49541 summarization stage: %s\n' "${FAILED}" >&2
  fi
}
trap cleanup EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite GSE49541 summarization artifact: %s\n' "${path}" >&2
    exit 1
  fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/validation" "${STAGE}/summary"
install -d -m 0700 "${CEL_WORK}"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}"           "${RUNTIME}/ARTIFACTS.json" \
  "${FIREWALL_SHA256}"          "${FIREWALL}/ARTIFACTS.json" \
  "${CROSSWALK_SHA256}"         "${CROSSWALK}/ARTIFACTS.json" \
  "${CROSSWALK_TSV_SHA256}"     "${CROSSWALK_TSV}" \
  "${RAW_AUDIT_SHA256}"         "${RAW_AUDIT}/ARTIFACTS.json" \
  "${CEL_MANIFEST_SHA256}"      "${CEL_MANIFEST}" \
  "${PLATFORM_AUDIT_SHA256}"    "${PLATFORM_AUDIT}/ARTIFACTS.json" \
  "${PLATFORM_AXIS_SHA256}"     "${PLATFORM_AXIS}" \
  "${PARTICIPANT_AUDIT_SHA256}" "${PARTICIPANT_AUDIT}/ARTIFACTS.json" \
  "${PARTICIPANTS_SHA256}"      "${PARTICIPANTS}" \
  "${PREPROCESSING_SHA256}"     "${PREPROCESSING}" \
  "${TASKSPEC_SHA256}"          "${TASKSPEC}" \
  "${GATE_SHA256}"              "${GATE}" \
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
  tests.unit.test_verify_gse49541_summarization_receipt \
  tests.unit.test_build_gse49541_participant_masks \
  tests.unit.test_stage_gse49541_cel_arrays \
  tests.unit.test_build_gpl570_gencode_v49_crosswalk \
  tests.contract.test_microarray_transfer_contract \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

## ---------------------------------- the separation verdict must already be a pass
python3 - "${FIREWALL}" <<'PY'
from pathlib import Path
import json
import sys

verdict = json.loads(
    (Path(sys.argv[1]) / "firewall_verdict.json").read_text(encoding="utf-8")
)
if (
    verdict["status"] != "pass_inductive_firewall_on_real_CEL"
    or verdict["series"] != "GSE49541"
    or verdict["method"] != "frma_with_exact_hgu133plus2frmavecs"
    or verdict["held_array_application"] != "one_array_at_a_time_with_frozen_training_object"
    or verdict["arrays_per_frma_call"] != 1
    or verdict["held_array_summary_invariant_to_other_arrays_in_the_run"] is not True
    or verdict["held_array_summary_invariant_to_processing_order"] is not True
    or verdict["perturbation_control_detected_for_every_target"] is not True
    or verdict["distinct_targets_have_distinct_digests"] is not True
    or verdict["multi_array_frma_call_is_batch_dependent"] is not True
    or verdict["multi_array_frma_call_used_as_preprocessing_path"] is not False
):
    raise SystemExit("refusing to summarize without a proven single-array firewall")
PY

## -------------------------------- stage the whole cohort onto node-local disk
python3 scripts/stage_gse49541_cel_arrays.py \
  --tar "${CEL_TAR}" \
  --manifest "${CEL_MANIFEST}" \
  --output "${CEL_WORK}" \
  --receipt "${STAGE}/summary/cel_staging_receipt.json" \
  > "${STAGE}/validation/cel_staging_stdout.json" \
  2> "${STAGE}/validation/cel_staging_stderr.txt"

## ------------------------------- one array at a time, never a multi-array batch
module purge
module load "${R_MODULE}"
R --version > "${STAGE}/validation/R.version.txt"
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/summarize_gse49541_gpl570_single_array.R" \
  "${CEL_WORK}" \
  "${CEL_MANIFEST}" \
  "${PLATFORM_AXIS}" \
  "${CROSSWALK_TSV}" \
  "${STAGE}/summary" \
  > "${STAGE}/validation/summarization_stdout.txt" \
  2> "${STAGE}/validation/summarization_stderr.txt"
cp "${STAGE}/summary/R_sessionInfo.txt" "${STAGE}/validation/R_sessionInfo.txt"
cp "${STAGE}/summary/package_versions.tsv" "${STAGE}/validation/package_versions.tsv"

## --------------------------------------------- masks and fail-closed adjudication
module purge
module load "${PYTHON_MODULE}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 scripts/build_gse49541_participant_masks.py \
  --participants "${PARTICIPANTS}" \
  --summarized-columns "${STAGE}/summary/summarized_columns.txt" \
  --output "${STAGE}/summary/participant_masks.json" \
  > "${STAGE}/validation/participant_masks_stdout.json" \
  2> "${STAGE}/validation/participant_masks_stderr.txt"

python3 scripts/verify_gse49541_summarization_receipt.py \
  --receipt "${STAGE}/summary/summarization_receipt.json" \
  --platform-matrix "${STAGE}/summary/gse49541_gpl570_platform_feature_matrix.tsv" \
  --gene-matrix "${STAGE}/summary/gse49541_gpl570_gene_matrix.tsv" \
  --output "${STAGE}/summary/summarization_verdict.json" \
  > "${STAGE}/validation/summarization_verdict_stdout.json" \
  2> "${STAGE}/validation/summarization_verdict_stderr.txt"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
rm -rf "${CEL_WORK}"

sha256sum \
  scripts/stage_gse49541_cel_arrays.py \
  scripts/summarize_gse49541_gpl570_single_array.R \
  scripts/build_gse49541_participant_masks.py \
  scripts/verify_gse49541_summarization_receipt.py \
  tests/unit/test_verify_gse49541_summarization_receipt.py \
  tests/unit/test_build_gse49541_participant_masks.py \
  slurm/summarize_gse49541_gpl570_single_array_cpu.sbatch \
  config/evaluation/microarray_transfer_preprocessing.toml \
  config/evaluation/gse49541_fibrosis_transfer_task.toml \
  config/evaluation/gse49541_fibrosis_transfer_promotion_gate.json \
  "${RUNTIME}/ARTIFACTS.json" \
  "${FIREWALL}/ARTIFACTS.json" \
  "${CROSSWALK}/ARTIFACTS.json" \
  "${RAW_AUDIT}/ARTIFACTS.json" \
  "${PARTICIPANT_AUDIT}/ARTIFACTS.json" \
  "${PLATFORM_AXIS}" \
  > "${STAGE}/source.sha256"

if find "${STAGE}" -type l -print -quit | grep -q .; then
  printf 'GSE49541 summarization stage contains a symlink\n' >&2
  exit 1
fi
if find "${STAGE}" -type f -name '*.CEL' -print -quit | grep -q .; then
  printf 'GSE49541 summarization stage contains a raw CEL file\n' >&2
  exit 1
fi

## ---------------------------------------------------------------------- freeze
python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1]).resolve(strict=True)
summary = stage / "summary"
verdict = json.loads((summary / "summarization_verdict.json").read_text(encoding="utf-8"))
masks = json.loads((summary / "participant_masks.json").read_text(encoding="utf-8"))
staging = json.loads((summary / "cel_staging_receipt.json").read_text(encoding="utf-8"))

if (
    verdict["status"] != "pass_label_blind_single_array_summarization"
    or verdict["arrays_summarized"] != 72
    or verdict["arrays_read_per_frma_call"] != 1
    or verdict["platform_features"] != 54_675
    or verdict["gene_features"] < 1
    or verdict["samples_excluded"] != 0
    or verdict["label_blind"] is not True
    or verdict["automatic_exclusion_from_descriptive_flag"] is not False
    or verdict["cross_array_QC_visible_to_training"] is not False
    or verdict["all_sample_RMA_run"] is not False
    or verdict["across_array_quantile_normalization_run"] is not False
    or verdict["GEO_series_matrix_read"] is not False
    or verdict["labels_read"] is not False
    or verdict["model_training_activated"] is not False
    or verdict["sealed_outcomes_read"] is not False
):
    raise SystemExit("GSE49541 summarization verdict differs from the frozen contract")
if (
    masks["participants"] != 72
    or masks["structurally_missing_imputed"] is not False
    or masks["outcome_columns_read"] is not False
    or masks["labels_read"] is not False
    or sorted(masks["structurally_missing_fields"]) != ["age", "exact_fibrosis_stage", "sex"]
):
    raise SystemExit("GSE49541 participant masks differ from the frozen contract")
if staging["staged_arrays"] != 72 or staging["every_array_matches_frozen_sha256"] is not True:
    raise SystemExit("GSE49541 CEL staging receipt differs")

metadata = {
    "artifact_class": "gse49541_gpl570_single_array_expression_matrices",
    "series": "GSE49541",
    "platform_id": "GPL570",
    "cohort_family_id": "gse31803_gse49541_fibrosis_array",
    "method": "frma_with_exact_hgu133plus2frmavecs",
    "arrays_summarized": verdict["arrays_summarized"],
    "arrays_read_per_frma_call": 1,
    "platform_features": verdict["platform_features"],
    "gene_features": verdict["gene_features"],
    "eligible_platform_features": verdict["eligible_platform_features"],
    "preprocessing_object_sha256": verdict["preprocessing_object_sha256"],
    "platform_feature_matrix_sha256": verdict["platform_feature_matrix_sha256"],
    "gene_matrix_sha256": verdict["gene_matrix_sha256"],
    "package_versions": verdict["package_versions"],
    "participants": masks["participants"],
    "structurally_missing_fields": masks["structurally_missing_fields"],
    "structurally_missing_imputed": False,
    "samples_excluded": 0,
    "automatic_exclusion_from_descriptive_flag": False,
    "cross_array_QC_visible_to_training": False,
    "label_blind": True,
    "all_sample_RMA_run": False,
    "across_array_quantile_normalization_run": False,
    "GEO_series_matrix_read": False,
    "raw_CEL_redistributed": False,
    "derived_matrix_redistribution": "source_and_model_specific_review_required",
    "labels_read": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
}
freeze_tree(summary, metadata)
verify_frozen_tree(summary)
freeze_tree(stage, metadata | {"artifact_class": "gse49541_gpl570_summarization_campaign"})
verify_frozen_tree(stage)
PY

trap - EXIT
rm -rf "${CEL_WORK}"
mv "${STAGE}" "${TARGET}"
printf 'GSE49541 GPL570 single-array summarization: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/summary/ARTIFACTS.json"
