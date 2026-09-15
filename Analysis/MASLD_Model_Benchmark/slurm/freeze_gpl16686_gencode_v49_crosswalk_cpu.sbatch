#!/usr/bin/env bash
#SBATCH --job-name=model-data-087
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Freeze the GPL16686 transcript cluster -> Entrez -> one-to-one Ensembl ->
# GENCODE v49 crosswalk, on the user's explicit mapping authority:
# hugene20sttranscriptcluster.db -> org.Hs.eg.db -> GENCODE v49.
#
# Two questions are settled here rather than assumed. The runtime probe reported
# 0 of 200 transcript clusters carrying an Entrez ID; this queries all of them.
# And the frozen requirements record platform_rows = 53981 from the GEO table while
# the core summarizer returns 53617, so the axis is taken from the annotation
# package, which is what the matrix will be keyed on.
#
# The Thermo Release 36 vendor asset stays recorded as unavailable provenance
# and is a model-card limitation. Coordinates and GB_ACC are never used.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export R_LIBS_USER=/nonexistent/masld_bench_r_user_library

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/microarray-r-bioc-21104442
RUNTIME_SHA256=8115e0966804dfc4b7bfac3f95016896d099d82fc51f3428fab1ec34182a7d07
FIREWALL=${ROOT}/executions/model-check-091-21109479-firewall/firewall
FIREWALL_SHA256=414b3d39bde46150fdf2b8205e29b6e85327ebb971f165551b95290455be055b
PREPROCESSING=${ROOT}/config/evaluation/microarray_transfer_preprocessing.toml
PREPROCESSING_SHA256=11162c381017bd16736433f478c0c82d97726ae2d896957d27efa5aee899c27e
LABEL_SEMANTICS=${ROOT}/config/evaluation/microarray_transfer_label_semantics.toml
LABEL_SEMANTICS_SHA256=b66c20ef5c6bc2ce7597affe7816358bd7fd4ca65df8e3f4be1096dbddaf0513
GENCODE_GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
GENCODE_GTF_SHA256=73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9

TARGET=${ROOT}/executions/model-data-087-${SLURM_JOB_ID}-crosswalk
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi }
trap preserve EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then printf 'Refusing to overwrite: %s\n' "${path}" >&2; exit 1; fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/validation" "${STAGE}/annotation"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}"         "${RUNTIME}/ARTIFACTS.json" \
  "${FIREWALL_SHA256}"        "${FIREWALL}/ARTIFACTS.json" \
  "${PREPROCESSING_SHA256}"   "${PREPROCESSING}" \
  "${LABEL_SEMANTICS_SHA256}" "${LABEL_SEMANTICS}" \
  "${GENCODE_GTF_SHA256}"     "${GENCODE_GTF}" \
  | sha256sum --check --strict

cd "${ROOT}"

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 -m unittest -v \
  tests.unit.test_build_gpl16686_gencode_v49_crosswalk \
  tests.contract.test_microarray_transfer_contract \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

module purge
module load R/4.5.2-mba
R --version > "${STAGE}/validation/R.version.txt"
R_LIBS="${RUNTIME}/library" Rscript --vanilla \
  "${ROOT}/scripts/dump_gpl16686_annotation_maps.R" "${STAGE}/annotation" \
  > "${STAGE}/validation/annotation_dump_stdout.txt" \
  2> "${STAGE}/validation/annotation_dump_stderr.txt"
cp "${STAGE}/annotation/R_sessionInfo.txt" "${STAGE}/validation/R_sessionInfo.txt"

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python3 scripts/build_gpl16686_gencode_v49_crosswalk.py \
  --platform-axis "${STAGE}/annotation/gpl16686_platform_axis.tsv" \
  --probe-entrez "${STAGE}/annotation/gpl16686_probe_entrez.tsv" \
  --entrez-ensembl "${STAGE}/annotation/entrez_ensembl.tsv" \
  --annotation-absent "${STAGE}/annotation/gpl16686_features_absent_from_annotation.tsv" \
  --entrez-absent "${STAGE}/annotation/entrez_absent_from_org_hs_eg_db.tsv" \
  --gencode-gtf "${GENCODE_GTF}" \
  --annotation-versions "${STAGE}/annotation/annotation_package_versions.tsv" \
  --output "${STAGE}/crosswalk" \
  > "${STAGE}/validation/crosswalk_stdout.json" \
  2> "${STAGE}/validation/crosswalk_stderr.txt"
cp "${STAGE}/annotation/annotation_package_versions.tsv" "${STAGE}/validation/package_versions.tsv"
python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"

sha256sum \
  scripts/dump_gpl16686_annotation_maps.R \
  scripts/build_gpl16686_gencode_v49_crosswalk.py \
  tests/unit/test_build_gpl16686_gencode_v49_crosswalk.py \
  slurm/freeze_gpl16686_gencode_v49_crosswalk_cpu.sbatch \
  config/evaluation/microarray_transfer_preprocessing.toml \
  config/evaluation/microarray_transfer_label_semantics.toml \
  "${RUNTIME}/ARTIFACTS.json" \
  "${FIREWALL}/ARTIFACTS.json" \
  "${GENCODE_GTF}" \
  > "${STAGE}/source.sha256"

python3 - "${STAGE}" "${GENCODE_GTF_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

stage = Path(sys.argv[1]).resolve(strict=True)
gencode_sha256 = sys.argv[2]
crosswalk = stage / "crosswalk"
receipt = json.loads((crosswalk / "crosswalk_receipt.json").read_text(encoding="utf-8"))
dump = json.loads((stage / "annotation/annotation_dump_receipt.json").read_text(encoding="utf-8"))

versions = receipt["package_versions"]
if (
    receipt["status"] != "pass_gpl16686_gencode_v49_crosswalk_frozen"
    or receipt["platform_id"] != "GPL16686"
    or receipt["platform_features"] != 53_617
    or receipt["gencode_gtf_sha256"] != gencode_sha256
    or receipt["annotation"] != "GENCODE_v49"
    or receipt["ambiguous_mapping_state"] != "join_unresolved"
    or receipt["retired_or_absent_GENCODE_state"] != "below_qc"
    or versions.get("hugene20sttranscriptcluster.db") != "8.8.0"
    or versions.get("org.Hs.eg.db") != "3.22.0"
    or versions.get("AnnotationDbi") != "1.72.0"
    or receipt["coordinate_or_GB_ACC_mapping_inference_used"] is not False
    or receipt["outcome_selected_probe_mapping"] is not False
    or receipt["depends_on_CEL_intensity"] is not False
    or receipt["depends_on_held_cohort_arrays"] is not False
    or receipt["depends_on_labels"] is not False
    or receipt["labels_read"] is not False
    or receipt["model_training_activated"] is not False
    or receipt["sealed_outcomes_read"] is not False
):
    raise SystemExit("GPL16686 crosswalk receipt differs from the frozen contract")
if dump["coordinate_or_GB_ACC_mapping_inference_used"] is not False:
    raise SystemExit("GPL16686 annotation dump inferred a mapping from coordinates")
if sum(receipt["contract_state_counts"].values()) != 53_617:
    raise SystemExit("GPL16686 crosswalk lost a platform feature")
if receipt["eligible_platform_features"] < 1 or receipt["one_to_one_genes"] < 1:
    raise SystemExit("GPL16686 crosswalk retained no one-to-one feature")

metadata = {
    "artifact_class": "gpl16686_gencode_v49_feature_crosswalk",
    "platform_id": "GPL16686",
    "series": "GSE83452",
    "cohort_family_id": "antwerp_inserm_shared",
    "annotation": "GENCODE_v49",
    "genome_build": "GRCh38.p14",
    "package_versions": versions,
    "mapping_authority": "hugene20sttranscriptcluster.db_then_org.Hs.eg.db_then_GENCODE_v49",
    "thermo_release_36_asset": "unavailable_permission_login_gated",
    "platform_features": receipt["platform_features"],
    "geo_platform_table_rows_not_used_as_axis": 53_981,
    "eligible_platform_features": receipt["eligible_platform_features"],
    "one_to_one_genes": receipt["one_to_one_genes"],
    "contract_state_counts": receipt["contract_state_counts"],
    "crosswalk_sha256": receipt["crosswalk_sha256"],
    "gencode_v49_crosswalk_frozen": True,
    "features_with_entrez": dump["features_with_entrez"],
    "first_200_keys_with_entrez": dump["first_200_keys_with_entrez"],
    "depends_on_CEL_intensity": False,
    "depends_on_held_cohort_arrays": False,
    "depends_on_labels": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed",
}
freeze_tree(crosswalk, metadata)
verify_frozen_tree(crosswalk)
freeze_tree(stage, metadata | {"artifact_class": "gpl16686_gencode_v49_crosswalk_campaign"})
verify_frozen_tree(stage)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GPL16686 GENCODE v49 crosswalk: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/crosswalk/ARTIFACTS.json"
