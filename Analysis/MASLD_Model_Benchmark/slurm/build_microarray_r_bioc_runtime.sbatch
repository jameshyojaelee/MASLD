#!/usr/bin/env bash
#SBATCH --job-name=model-runtime-145
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Frozen portable Bioconductor 3.22 runtime for the two blocked microarray
# transfer TaskSpecs (GSE49541 / GPL570 and GSE83452 / GPL16686).
#
# The runtime is an installed R library overlay, not a conda prefix, because a
# conda prefix hardcodes R_HOME_DIR into bin/R and therefore cannot survive the
# staging -> target move that the freeze protocol requires. An installed R
# library tree is relocatable. R itself is the fixed R/4.5.2-mba module, which
# already ships Bioconductor 3.22 AnnotationDbi 1.72.0, Biostrings 2.78.0 and
# BiocManager 1.30.27. This matches the seurat-5.5.1-portable runtime shape.
#
# This job installs software and runs a capability probe on ONE array per
# platform. It never runs an all-sample RMA, never runs an across-array quantile
# normalization, never reads a GEO series matrix, never reads a label, and never
# fits a model.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER=/dev/null
export R_LIBS_USER=/nonexistent/masld_bench_r_user_library

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PREPARATION=${ROOT}/executions/model-data-072-21080651-preparation
PREPARATION_ARTIFACTS_SHA256=8dd6f7288f7c06d7e1db5e64d87eedb3cfa122c13afee3da6293ff204144f96b
PREPROCESSING=${ROOT}/config/evaluation/microarray_transfer_preprocessing.toml
PREPROCESSING_SHA256=11162c381017bd16736433f478c0c82d97726ae2d896957d27efa5aee899c27e
R_MODULE=R/4.5.2-mba
BIOC_RELEASE=3.22

TARGET=${ROOT}/executions/environments/microarray-r-bioc-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed microarray runtime stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
for path in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite microarray Bioconductor runtime: %s\n' "${path}" >&2
    exit 1
  fi
done

install -d -m 0750 "${STAGE}" "${STAGE}/library" "${STAGE}/downloads" "${STAGE}/validation"
printf '%s  %s\n' \
  "${PREPARATION_ARTIFACTS_SHA256}" "${PREPARATION}/ARTIFACTS.json" \
  "${PREPROCESSING_SHA256}" "${PREPROCESSING}" \
  | sha256sum --check --strict

GPL570_FIXTURE=$(find "${PREPARATION}/fixtures" -maxdepth 1 -type f -name 'GSE49541_*.CEL' -print -quit)
GPL16686_FIXTURE=$(find "${PREPARATION}/fixtures" -maxdepth 1 -type f -name 'GSE83452_*.CEL' -print -quit)
if [[ -z "${GPL570_FIXTURE}" || -z "${GPL16686_FIXTURE}" ]]; then
  printf 'Frozen single-array CEL fixture is absent.\n' >&2
  exit 1
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

## --------------------------------------------- portable Bioconductor overlay
module purge
module load "${R_MODULE}"
R --version > "${STAGE}/validation/R.version.txt"
Rscript --vanilla "${ROOT}/scripts/install_microarray_bioc_library.R" \
  "${STAGE}/library" \
  "${STAGE}/downloads" \
  "${STAGE}/validation/package_manifest.tsv" \
  > "${STAGE}/validation/install.stdout.txt" \
  2> "${STAGE}/validation/install.stderr.txt"

## ---------------------------------------------------- single-array capability
R_LIBS="${STAGE}/library" Rscript --vanilla \
  "${ROOT}/scripts/probe_microarray_single_array_runtime.R" \
  "${GPL570_FIXTURE}" "${GPL16686_FIXTURE}" "${STAGE}/validation"

## ------------------------------ 55-key runtime receipt for the audit requirements
R_LIBS="${STAGE}/library" Rscript --vanilla \
  "${ROOT}/scripts/audit_microarray_r_runtime.R" \
  microarray_r_bioc_portable \
  "${GPL570_FIXTURE}" "${GPL16686_FIXTURE}" \
  "${STAGE}/validation/microarray_r_bioc_portable.tsv"

## ------------------------------------------------------------------ provenance
sha256sum "${STAGE}"/downloads/*.tar.gz > "${STAGE}/validation/source-tarballs.sha256"
cd "${ROOT}"
sha256sum \
  slurm/build_microarray_r_bioc_runtime.sbatch \
  scripts/install_microarray_bioc_library.R \
  scripts/probe_microarray_single_array_runtime.R \
  scripts/audit_microarray_r_runtime.R \
  "${PREPROCESSING}" \
  "${PREPARATION}/ARTIFACTS.json" \
  "${GPL570_FIXTURE}" \
  "${GPL16686_FIXTURE}" \
  > "${STAGE}/source.sha256"

if find "${STAGE}" -type l -print -quit | grep -q .; then
  printf 'Microarray Bioconductor runtime contains a symlink\n' >&2
  find "${STAGE}" -type l -print >&2
  exit 1
fi

## ---------------------------------------------------------------------- freeze
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${STAGE}" "${BIOC_RELEASE}" "${R_MODULE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads(
    (root / "validation/capability_probe.json").read_text(encoding="utf-8")
)
capabilities = probe["capabilities"]
if (
    probe.get("status") != "pass_single_array_runtime_capability"
    or capabilities.get("calvin_header_parse_both_platforms") is not True
    or capabilities.get("GPL570_single_array_summarization") is not True
    or capabilities.get("GPL16686_single_array_summarization") is not True
    or capabilities.get("probe_to_entrez_to_ensembl_both_platforms") is not True
    or probe.get("all_sample_RMA_run") is not False
    or probe.get("across_array_quantile_normalization_run") is not False
    or probe.get("GEO_series_matrix_read") is not False
    or probe.get("expression_values_exported") is not False
    or probe.get("labels_read") is not False
    or probe.get("model_training_activated") is not False
    or probe.get("sealed_outcomes_read") is not False
):
    raise SystemExit("microarray single-array capability probe differs")

required = {
    "affy": "1.88.0",
    "affyio": "1.80.0",
    "affxparser": "1.82.0",
    "oligo": "1.74.0",
    "frma": "1.62.0",
    "SCAN.UPC": "2.52.0",
    "AnnotationDbi": "1.72.0",
    "org.Hs.eg.db": "3.22.0",
    "hgu133plus2cdf": "2.18.0",
    "hgu133plus2.db": "3.13.0",
    "hgu133plus2frmavecs": "1.5.0",
    "hugene20sttranscriptcluster.db": "8.8.0",
    "pd.hugene.2.0.st": "3.14.1",
    "pd.hg.u133.plus.2": "3.12.0",
}
observed = {name: probe["package_versions"].get(name) for name in required}
if observed != required:
    raise SystemExit(f"microarray runtime package versions differ: {observed!r}")

with (root / "validation/microarray_r_bioc_portable.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    receipt = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
for key in (
    "capability.Calvin_header_reader_both_platforms",
    "capability.GPL570_exact_frma",
    "capability.GPL570_SCAN",
    "capability.GPL570_affy_summary",
    "capability.GPL16686_SCAN",
    "capability.GPL16686_oligo_core_summary",
    "capability.GPL16686_transcript_cluster_annotation",
    "capability.GPL570_expression_rank_baseline",
    "capability.GPL16686_expression_rank_baseline",
):
    if receipt.get(key) != "true":
        raise SystemExit(f"microarray runtime receipt capability differs: {key}")

freeze_tree(
    root,
    {
        "artifact_class": "microarray_gpl570_gpl16686_single_array_bioconductor_runtime",
        "runtime_id": "microarray_r_bioc_portable",
        "runtime_shape": "relocatable_r_library_overlay_on_fixed_r_module",
        "r_module": sys.argv[3],
        "R_version": probe["R_version"],
        "bioconductor_release": sys.argv[2],
        "package_versions": required,
        "GPL570_single_array_method": "frma_with_exact_hgu133plus2frmavecs",
        "GPL16686_single_array_method": "SCAN_or_normalize_free_core_summary",
        "GPL16686_vendor_release_36_asset": "unavailable_permission_login_gated",
        "GPL16686_mapping_authority": (
            "hugene20sttranscriptcluster.db_then_org.Hs.eg.db_then_GENCODE_v49"
        ),
        "arrays_read_per_platform": 1,
        "dependency_install_performed": True,
        "all_sample_RMA_run": False,
        "across_array_quantile_normalization_run": False,
        "GEO_series_matrix_read": False,
        "expression_values_exported": False,
        "gencode_v49_crosswalk_frozen": False,
        "labels_read": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Microarray Bioconductor single-array runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
