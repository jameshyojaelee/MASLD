#!/usr/bin/env bash
#SBATCH --job-name=model-data-088
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Merge the eight GSE83452 shards into the cohort matrices.
#
# The merge is only legitimate if an array's summary does not depend on which
# job produced it, so the sentinel arrays that appear in every shard must agree
# bitwise across all eight independent jobs before anything is assembled. Each
# column is also re-digested from the shard file and matched against the digest
# that shard recorded, so a column cannot be transposed or dropped silently.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PLAN=${ROOT}/executions/gse83452-shard-plan-20260825/shard_plan.json
CROSSWALK=${ROOT}/executions/model-data-087-21111875-crosswalk/crosswalk
CROSSWALK_SHA256=c6bfe6c3d44b98fd612cf81ba48014fb32968a4d690b9424369c7bc23018048d
FIREWALL=${ROOT}/executions/model-check-091-21109479-firewall/firewall
FIREWALL_SHA256=414b3d39bde46150fdf2b8205e29b6e85327ebb971f165551b95290455be055b
# microarray_transfer_preprocessing.toml pins model-data-067 as
# participant_audit_v2. That table predates the cross-accession overlap audit
# (model-data-072) and still records cohort_family_id = gse83452_longitudinal_array,
# which would UNDER-group the family and let GSE106737 be treated as an
# independent cohort. The corrected table from model-data-072 is used instead;
# the stale requirements pin is reported, not silently reconciled.
RECORDS=${ROOT}/executions/model-data-072-21080651-participants/gse83452_records.tsv
RECORDS_SHA256=04cdbdfa018c379eb708a93a73de419cb848da81a3a26b33560c204b6bacac94
PARTICIPANT_AUDIT=${ROOT}/executions/model-data-072-21080651-participants
PARTICIPANT_AUDIT_SHA256=6655f23b2092e163532113e74d75d3c2a93857ee870b4f0a94fdcb4fa7f9dbc0
OVERLAP_AUDIT=${ROOT}/executions/model-data-072-21080651-overlap
OVERLAP_AUDIT_SHA256=e577e1f581a484095748a888cfaebcfb83548d452bab523eab8ef19dbb0c89eb
LABEL_SEMANTICS=${ROOT}/config/evaluation/microarray_transfer_label_semantics.toml
LABEL_SEMANTICS_SHA256=b66c20ef5c6bc2ce7597affe7816358bd7fd4ca65df8e3f4be1096dbddaf0513
BASELINE_TASK=${ROOT}/config/evaluation/gse83452_baseline_nash_transfer_task.toml
BASELINE_TASK_SHA256=c567022d32d478b364d6a3cd9c84ba692a023bac930d99d3bb262499d60e943b
PAIRED_TASK=${ROOT}/config/evaluation/gse83452_paired_intervention_stress_task.toml
PAIRED_TASK_SHA256=da8e5e8f5d58767905f2f3862d6cd22bb80f25d1bfacdb17f28f9fd2586b4d75

TARGET=${ROOT}/executions/model-data-088-${SLURM_JOB_ID}-summarization
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi }
trap preserve EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then printf 'Refusing to overwrite: %s\n' "${path}" >&2; exit 1; fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/validation"

printf '%s  %s\n' \
  "${CROSSWALK_SHA256}"         "${CROSSWALK}/ARTIFACTS.json" \
  "${FIREWALL_SHA256}"          "${FIREWALL}/ARTIFACTS.json" \
  "${RECORDS_SHA256}"           "${RECORDS}" \
  "${PARTICIPANT_AUDIT_SHA256}" "${PARTICIPANT_AUDIT}/ARTIFACTS.json" \
  "${OVERLAP_AUDIT_SHA256}"     "${OVERLAP_AUDIT}/ARTIFACTS.json" \
  "${LABEL_SEMANTICS_SHA256}"   "${LABEL_SEMANTICS}" \
  "${BASELINE_TASK_SHA256}"     "${BASELINE_TASK}" \
  "${PAIRED_TASK_SHA256}"       "${PAIRED_TASK}" \
  | sha256sum --check --strict

cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"

python3 -m unittest -v \
  tests.unit.test_gse83452_shard_plan_and_merge \
  tests.unit.test_build_gse83452_participant_masks \
  tests.unit.test_build_gpl16686_gencode_v49_crosswalk \
  tests.unit.test_gse83452_oligo_firewall_receipt \
  tests.contract.test_microarray_transfer_contract \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

SHARD_ARGS=()
for d in "${ROOT}"/executions/data-stage-260-*-shard*/shard; do
  SHARD_ARGS+=(--shard-dir "${d}")
done
if [[ ${#SHARD_ARGS[@]} -ne 16 ]]; then
  printf 'Expected 8 shard directories, found %d\n' $((${#SHARD_ARGS[@]} / 2)) >&2
  exit 1
fi

python3 scripts/merge_gse83452_gpl16686_shards.py \
  "${SHARD_ARGS[@]}" --plan "${PLAN}" --output "${STAGE}/summary" \
  > "${STAGE}/validation/merge_stdout.json" 2> "${STAGE}/validation/merge_stderr.txt"

python3 scripts/build_gse83452_participant_masks.py \
  --records "${RECORDS}" \
  --summarized-columns "${STAGE}/summary/summarized_columns.txt" \
  --output "${STAGE}/summary/participant_masks.json" \
  > "${STAGE}/validation/participant_masks_stdout.json" \
  2> "${STAGE}/validation/participant_masks_stderr.txt"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
cp "${CROSSWALK}/crosswalk_receipt.json" "${STAGE}/validation/crosswalk_receipt.json"

sha256sum \
  scripts/merge_gse83452_gpl16686_shards.py \
  scripts/build_gse83452_participant_masks.py \
  scripts/plan_gse83452_summarization_shards.py \
  scripts/summarize_gse83452_gpl16686_shard.R \
  tests/unit/test_gse83452_shard_plan_and_merge.py \
  tests/unit/test_build_gse83452_participant_masks.py \
  slurm/merge_gse83452_gpl16686_shards_cpu.sbatch \
  config/evaluation/microarray_transfer_label_semantics.toml \
  config/evaluation/gse83452_baseline_nash_transfer_task.toml \
  config/evaluation/gse83452_paired_intervention_stress_task.toml \
  "${CROSSWALK}/ARTIFACTS.json" "${FIREWALL}/ARTIFACTS.json" "${RECORDS}" \
  "${OVERLAP_AUDIT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1]).resolve(strict=True)
summary = stage / "summary"
receipt = json.loads((summary / "merged_receipt.json").read_text(encoding="utf-8"))
masks = json.loads((summary / "participant_masks.json").read_text(encoding="utf-8"))

if (
    receipt["status"] != "pass_label_blind_single_array_summarization"
    or receipt["arrays_summarized"] != 231
    or receipt["arrays_read_per_summarization_call"] != 1
    or receipt["platform_features"] != 53_617
    or receipt["sentinel_digests_agree_across_shards"] is not True
    or receipt["cross_job_single_array_reproducibility_verified"] is not True
    or receipt["every_column_digest_matched_its_shard_receipt"] is not True
    or receipt["matrices_round_trip_bitwise"] is not True
    or receipt["samples_excluded"] != 0
    or receipt["label_blind"] is not True
    or receipt["all_sample_RMA_run"] is not False
    or receipt["across_array_quantile_normalization_run"] is not False
    or receipt["GEO_series_matrix_read"] is not False
    or receipt["labels_read"] is not False
    or receipt["model_training_activated"] is not False
    or receipt["sealed_outcomes_read"] is not False
):
    raise SystemExit("GSE83452 merged receipt differs from the frozen contract")
if (
    masks["records"] != 231
    or masks["participants"] != 171
    or masks["withheld_fields_imputed"] is not False
    or masks["outcome_columns_read"] is not False
    or masks["labels_read"] is not False
    or masks["independent_evaluation_between_accession_aliases"] is not False
    or masks["distinct_outer_groups"] != 171
):
    raise SystemExit("GSE83452 participant masks differ from the frozen contract")

metadata = {
    "artifact_class": "gse83452_gpl16686_single_array_expression_matrices",
    "series": "GSE83452",
    "platform_id": "GPL16686",
    "cohort_family_id": "antwerp_inserm_shared",
    "accession_aliases": ["GSE106737", "GSE83452"],
    "independent_evaluation_between_accession_aliases": False,
    "method": "single_array_core_summary_normalize_false",
    "arrays_summarized": receipt["arrays_summarized"],
    "arrays_read_per_summarization_call": 1,
    "platform_features": receipt["platform_features"],
    "gene_features": receipt["gene_features"],
    "shards": receipt["shards"],
    "sentinels": receipt["sentinels"],
    "cross_job_single_array_reproducibility_verified": True,
    "platform_feature_matrix_sha256": receipt["platform_feature_matrix_sha256"],
    "gene_matrix_sha256": receipt["gene_matrix_sha256"],
    "matrices_round_trip_bitwise": True,
    "participants": masks["participants"],
    "distinct_outer_groups": masks["distinct_outer_groups"],
    "outer_group_key_rule": masks["outer_group_key_rule"],
    "samples_excluded": 0,
    "label_blind": True,
    "all_sample_RMA_run": False,
    "across_array_quantile_normalization_run": False,
    "GEO_series_matrix_read": False,
    "raw_CEL_redistributed": False,
    "derived_matrix_redistribution": "source_and_model_specific_review_required",
    "thermo_release_36_asset": "unavailable_permission_login_gated",
    "labels_read": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
}
freeze_tree(summary, metadata)
verify_frozen_tree(summary)
freeze_tree(stage, metadata | {"artifact_class": "gse83452_gpl16686_summarization_campaign"})
verify_frozen_tree(stage)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE83452 GPL16686 single-array summarization: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/summary/ARTIFACTS.json"
