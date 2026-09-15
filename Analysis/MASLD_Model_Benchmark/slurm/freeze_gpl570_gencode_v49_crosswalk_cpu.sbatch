#!/usr/bin/env bash
#SBATCH --job-name=model-data-083
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Freeze the GPL570 platform feature -> Entrez -> one-to-one Ensembl ->
# GENCODE v49 crosswalk declared by microarray_transfer_preprocessing.toml.
#
# This closes the last declared gap on the GSE49541 lane,
# `gencode_v49_crosswalk_frozen: false`.
#
# The crosswalk is an annotation object. It is derived from the frozen GEO
# platform axis, the exact pinned annotation package versions, and GENCODE v49.
# It never opens a CEL file, never reads a GEO series matrix, never reads a
# label, and never fits a model, so it is identical whichever arrays are later
# held out. Ambiguity keeps the requirements' `join_unresolved` state and a gene
# GENCODE v49 retired keeps `below_qc`; neither is broken by a heuristic and
# neither is inferred from coordinates or GB_ACC.

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

PLATFORM_AUDIT=${ROOT}/executions/model-data-067-21080193-platform
PLATFORM_AUDIT_SHA256=0ce6af4196d2fb0905d6ca3d07fbbeb3bf1ce5b37776a63a726cc1d21794f4bf
PLATFORM_AXIS=${PLATFORM_AUDIT}/GPL570_selected_annotation.tsv
PLATFORM_AXIS_SHA256=4bbd91765d1a80b12439e8c39b2aca9218a4135a80fc6a50f6f4eed276d61b69

PREPROCESSING=${ROOT}/config/evaluation/microarray_transfer_preprocessing.toml
PREPROCESSING_SHA256=11162c381017bd16736433f478c0c82d97726ae2d896957d27efa5aee899c27e
TASKSPEC=${ROOT}/config/evaluation/gse49541_fibrosis_transfer_task.toml
TASKSPEC_SHA256=9b11584b3c1b2a01d4bb262790274032c31f92ff0b1da02a66588aa3e75a6ec3
GATE=${ROOT}/config/evaluation/gse49541_fibrosis_transfer_promotion_gate.json
GATE_SHA256=791dfa3199c3bdbdeb74bfae429cc8820364862c7187f5d1dbe11b46c03dc384

GENCODE_GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
GENCODE_GTF_SHA256=73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9

TARGET=${ROOT}/executions/model-data-083-${SLURM_JOB_ID}-crosswalk
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed GPL570 crosswalk stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite GPL570 crosswalk artifact: %s\n' "${path}" >&2
    exit 1
  fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/validation" "${STAGE}/annotation"

## ------------------------------------------- every frozen input is bound first
printf '%s  %s\n' \
  "${RUNTIME_SHA256}"         "${RUNTIME}/ARTIFACTS.json" \
  "${PLATFORM_AUDIT_SHA256}"  "${PLATFORM_AUDIT}/ARTIFACTS.json" \
  "${PLATFORM_AXIS_SHA256}"   "${PLATFORM_AXIS}" \
  "${PREPROCESSING_SHA256}"   "${PREPROCESSING}" \
  "${TASKSPEC_SHA256}"        "${TASKSPEC}" \
  "${GATE_SHA256}"            "${GATE}" \
  "${GENCODE_GTF_SHA256}"     "${GENCODE_GTF}" \
  | sha256sum --check --strict

cd "${ROOT}"

## ------------------------------------------------------------- unit tests
module purge
module load "${PYTHON_MODULE}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 -m unittest -v \
  tests.unit.test_build_gpl570_gencode_v49_crosswalk \
  tests.contract.test_microarray_transfer_contract \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

## ------------------------------- annotation authority dump (Bioconductor half)
module purge
module load "${R_MODULE}"
R --version > "${STAGE}/validation/R.version.txt"
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/dump_gpl570_annotation_maps.R" \
  "${PLATFORM_AXIS}" \
  "${STAGE}/annotation" \
  > "${STAGE}/validation/annotation_dump_stdout.txt" \
  2> "${STAGE}/validation/annotation_dump_stderr.txt"
cp "${STAGE}/annotation/R_sessionInfo.txt" "${STAGE}/validation/R_sessionInfo.txt"

## --------------------------------- GENCODE v49 intersection and state machine
module purge
module load "${PYTHON_MODULE}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 scripts/build_gpl570_gencode_v49_crosswalk.py \
  --platform-axis "${PLATFORM_AXIS}" \
  --probe-entrez "${STAGE}/annotation/gpl570_probe_entrez.tsv" \
  --entrez-ensembl "${STAGE}/annotation/entrez_ensembl.tsv" \
  --annotation-absent "${STAGE}/annotation/gpl570_features_absent_from_annotation.tsv" \
  --entrez-absent "${STAGE}/annotation/entrez_absent_from_org_hs_eg_db.tsv" \
  --gencode-gtf "${GENCODE_GTF}" \
  --annotation-versions "${STAGE}/annotation/annotation_package_versions.tsv" \
  --output "${STAGE}/crosswalk" \
  > "${STAGE}/validation/crosswalk_stdout.json" \
  2> "${STAGE}/validation/crosswalk_stderr.txt"
cp "${STAGE}/annotation/annotation_package_versions.tsv" \
  "${STAGE}/validation/package_versions.tsv"

python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"

sha256sum \
  scripts/dump_gpl570_annotation_maps.R \
  scripts/build_gpl570_gencode_v49_crosswalk.py \
  tests/unit/test_build_gpl570_gencode_v49_crosswalk.py \
  slurm/freeze_gpl570_gencode_v49_crosswalk_cpu.sbatch \
  config/evaluation/microarray_transfer_preprocessing.toml \
  config/evaluation/gse49541_fibrosis_transfer_task.toml \
  config/evaluation/gse49541_fibrosis_transfer_promotion_gate.json \
  "${RUNTIME}/ARTIFACTS.json" \
  "${PLATFORM_AUDIT}/ARTIFACTS.json" \
  "${PLATFORM_AXIS}" \
  "${GENCODE_GTF}" \
  > "${STAGE}/source.sha256"

if find "${STAGE}" -type l -print -quit | grep -q .; then
  printf 'GPL570 crosswalk stage contains a symlink\n' >&2
  find "${STAGE}" -type l -print >&2
  exit 1
fi

## ---------------------------------------------------------------------- freeze
python3 - "${STAGE}" "${GENCODE_GTF_SHA256}" "${PLATFORM_AXIS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1]).resolve(strict=True)
gencode_sha256 = sys.argv[2]
axis_sha256 = sys.argv[3]

crosswalk = stage / "crosswalk"
receipt = json.loads((crosswalk / "crosswalk_receipt.json").read_text(encoding="utf-8"))
dump = json.loads(
    (stage / "annotation/annotation_dump_receipt.json").read_text(encoding="utf-8")
)

versions = receipt["package_versions"]
if (
    receipt["status"] != "pass_gpl570_gencode_v49_crosswalk_frozen"
    or receipt["platform_id"] != "GPL570"
    or receipt["platform_features"] != 54_675
    or receipt["gencode_gtf_sha256"] != gencode_sha256
    or receipt["platform_axis_sha256"] != axis_sha256
    or receipt["annotation"] != "GENCODE_v49"
    or receipt["genome_build"] != "GRCh38.p14"
    or receipt["ambiguous_mapping_state"] != "join_unresolved"
    or receipt["retired_or_absent_GENCODE_state"] != "below_qc"
    or versions.get("hgu133plus2.db") != "3.13.0"
    or versions.get("org.Hs.eg.db") != "3.22.0"
    or versions.get("AnnotationDbi") != "1.72.0"
    or receipt["coordinate_or_GB_ACC_mapping_inference_used"] is not False
    or receipt["outcome_selected_probe_mapping"] is not False
    or receipt["depends_on_CEL_intensity"] is not False
    or receipt["depends_on_held_cohort_arrays"] is not False
    or receipt["depends_on_labels"] is not False
    or receipt["CEL_expression_values_read"] is not False
    or receipt["GEO_series_matrix_read"] is not False
    or receipt["labels_read"] is not False
    or receipt["model_training_activated"] is not False
    or receipt["sealed_outcomes_read"] is not False
):
    raise SystemExit("GPL570 crosswalk receipt differs from the frozen contract")
if (
    dump["ambiguity_resolved_in_this_script"] is not False
    or dump["CEL_expression_values_read"] is not False
    or dump["labels_read"] is not False
    or dump["platform_features"] != 54_675
):
    raise SystemExit("GPL570 annotation dump receipt differs")
if sum(receipt["contract_state_counts"].values()) != 54_675:
    raise SystemExit("GPL570 crosswalk lost a platform feature")
if receipt["eligible_platform_features"] < 1 or receipt["one_to_one_genes"] < 1:
    raise SystemExit("GPL570 crosswalk retained no one-to-one feature")

freeze_tree(
    crosswalk,
    {
        "artifact_class": "gpl570_gencode_v49_feature_crosswalk",
        "platform_id": "GPL570",
        "series": "GSE49541",
        "cohort_family_id": "gse31803_gse49541_fibrosis_array",
        "annotation": "GENCODE_v49",
        "genome_build": "GRCh38.p14",
        "package_versions": versions,
        "platform_features": receipt["platform_features"],
        "eligible_platform_features": receipt["eligible_platform_features"],
        "one_to_one_genes": receipt["one_to_one_genes"],
        "contract_state_counts": receipt["contract_state_counts"],
        "crosswalk_sha256": receipt["crosswalk_sha256"],
        "gencode_v49_crosswalk_frozen": True,
        "depends_on_CEL_intensity": False,
        "depends_on_held_cohort_arrays": False,
        "depends_on_labels": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(crosswalk)

freeze_tree(
    stage,
    {
        "artifact_class": "gpl570_gencode_v49_crosswalk_campaign",
        "platform_id": "GPL570",
        "series": "GSE49541",
        "gencode_v49_crosswalk_frozen": True,
        "crosswalk_sha256": receipt["crosswalk_sha256"],
        "eligible_platform_features": receipt["eligible_platform_features"],
        "one_to_one_genes": receipt["one_to_one_genes"],
        "model_training_activated": False,
        "labels_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(stage)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GPL570 GENCODE v49 crosswalk: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/crosswalk/ARTIFACTS.json"
