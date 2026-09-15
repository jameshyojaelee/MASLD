#!/usr/bin/env bash
#SBATCH --job-name=model-data-624
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Re-run the microarray runtime audit with the activated Bioconductor runtime in
# the inventory. Nothing about the three pre-existing runtimes changes; the new
# runtime is audited beside them with the same fixed 55-key receipt schema, and
# its receipt is produced from the frozen runtime's final path so that the audit
# also proves the library overlay survived the staging move.
#
# This job reads no label, fits no model, runs no all-sample RMA, and reads no
# GEO series matrix.

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
INDUCTIVE=${ROOT}/executions/model-data-072-21080651-inductive
INDUCTIVE_ARTIFACTS_SHA256=ea27d2ae95b914b1e01831798a37db692f681d504f6583606f0c902439a0cd58
PREPROCESSING=${ROOT}/config/evaluation/microarray_transfer_preprocessing.toml
PREPROCESSING_SHA256=11162c381017bd16736433f478c0c82d97726ae2d896957d27efa5aee899c27e
LABELS=${ROOT}/config/evaluation/microarray_transfer_label_semantics.toml
LABELS_SHA256=b66c20ef5c6bc2ce7597affe7816358bd7fd4ca65df8e3f4be1096dbddaf0513
RUNTIME_ENVIRONMENT=${ROOT}/executions/environments/microarray-r-bioc-__RUNTIME_JOB_ID__
RUNTIME_ARTIFACTS_SHA256=__RUNTIME_ARTIFACTS_SHA256__

RUN_STEM=${SLURM_JOB_NAME}-${SLURM_JOB_ID}
RUNTIME=${ROOT}/executions/${RUN_STEM}-runtime
RUNTIME_STAGE=${RUNTIME}.staging
RECEIPTS_STAGE=${RUNTIME}.receipts.staging
preserve_failed_stages() {
  for stage in "${RUNTIME_STAGE}" "${RECEIPTS_STAGE}"; do
    if [[ -e "${stage}" && ! -e "${stage}.failed" ]]; then
      mv "${stage}" "${stage}.failed"
    fi
  done
}
trap preserve_failed_stages EXIT
for path in "${RUNTIME}" "${RUNTIME_STAGE}" "${RECEIPTS_STAGE}"; do
  if [[ -e "${path}" || -e "${path}.failed" ]]; then
    printf 'Refusing to overwrite runtime audit artifact: %s\n' "${path}" >&2
    exit 1
  fi
done

cd "${ROOT}"
printf '%s  %s\n' \
  "${PREPARATION_ARTIFACTS_SHA256}" "${PREPARATION}/ARTIFACTS.json" \
  "${INDUCTIVE_ARTIFACTS_SHA256}" "${INDUCTIVE}/ARTIFACTS.json" \
  "${PREPROCESSING_SHA256}" "${PREPROCESSING}" \
  "${LABELS_SHA256}" "${LABELS}" \
  "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME_ENVIRONMENT}/ARTIFACTS.json" \
  | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python - "${PREPARATION}" "${INDUCTIVE}" "${RUNTIME_ENVIRONMENT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
python -m py_compile scripts/finalize_microarray_runtime_audit.py
python -m unittest -v tests.unit.test_finalize_microarray_runtime_audit

mkdir "${RECEIPTS_STAGE}"
GPL570_FIXTURE=$(find "${PREPARATION}/fixtures" -maxdepth 1 -type f -name 'GSE49541_*.CEL' -print -quit)
GPL16686_FIXTURE=$(find "${PREPARATION}/fixtures" -maxdepth 1 -type f -name 'GSE83452_*.CEL' -print -quit)
module purge
module load R/4.4.3
Rscript scripts/audit_microarray_r_runtime.R \
  module_R_4.4.3 "${GPL570_FIXTURE}" "${GPL16686_FIXTURE}" \
  "${RECEIPTS_STAGE}/module_R_4.4.3.tsv"
module purge
module load R/4.4.1
Rscript scripts/audit_microarray_r_runtime.R \
  module_R_4.4.1 "${GPL570_FIXTURE}" "${GPL16686_FIXTURE}" \
  "${RECEIPTS_STAGE}/module_R_4.4.1.tsv"
module purge
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript \
  scripts/audit_microarray_r_runtime.R \
  micromamba_rnaseq "${GPL570_FIXTURE}" "${GPL16686_FIXTURE}" \
  "${RECEIPTS_STAGE}/micromamba_rnaseq.tsv"
module load R/4.5.2-mba
R_LIBS="${RUNTIME_ENVIRONMENT}/library" Rscript --vanilla \
  scripts/audit_microarray_r_runtime.R \
  microarray_r_bioc_portable "${GPL570_FIXTURE}" "${GPL16686_FIXTURE}" \
  "${RECEIPTS_STAGE}/microarray_r_bioc_portable.tsv"

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python scripts/finalize_microarray_runtime_audit.py \
  --preparation "${PREPARATION}" \
  --inductive-fixture "${INDUCTIVE}" \
  --runtime-receipt "${RECEIPTS_STAGE}/module_R_4.4.3.tsv" \
  --runtime-receipt "${RECEIPTS_STAGE}/module_R_4.4.1.tsv" \
  --runtime-receipt "${RECEIPTS_STAGE}/micromamba_rnaseq.tsv" \
  --runtime-receipt "${RECEIPTS_STAGE}/microarray_r_bioc_portable.tsv" \
  --single-array-capability "${RUNTIME_ENVIRONMENT}/validation/capability_probe.json" \
  --dependency-install-performed \
  --preprocessing-contract "${PREPROCESSING}" \
  --label-contract "${LABELS}" \
  --output "${RUNTIME_STAGE}"
mv "${RECEIPTS_STAGE}"/*.tsv "${RUNTIME_STAGE}/"
rmdir "${RECEIPTS_STAGE}"
cp -p "${RUNTIME_ENVIRONMENT}/validation/capability_probe.json" \
  "${RUNTIME_STAGE}/single_array_capability_probe.json"
cp -p "${RUNTIME_ENVIRONMENT}/validation/package_manifest.tsv" \
  "${RUNTIME_STAGE}/runtime_package_manifest.tsv"
sha256sum \
  scripts/audit_microarray_r_runtime.R \
  scripts/finalize_microarray_runtime_audit.py \
  scripts/install_microarray_bioc_library.R \
  scripts/probe_microarray_single_array_runtime.R \
  tests/unit/test_finalize_microarray_runtime_audit.py \
  slurm/audit_microarray_runtime_activation.sbatch \
  "${PREPROCESSING}" "${LABELS}" \
  "${PREPARATION}/ARTIFACTS.json" "${INDUCTIVE}/ARTIFACTS.json" \
  "${RUNTIME_ENVIRONMENT}/ARTIFACTS.json" \
  > "${RUNTIME_STAGE}/source.sha256"
python - "${RUNTIME_STAGE}" "${RUNTIME_ENVIRONMENT}" "${RUNTIME_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "microarray_runtime_audit.json").read_text(encoding="utf-8"))
if (
    audit.get("synthetic_training_only_quantile_fixture_passed") is not True
    or audit.get("synthetic_fixture_is_biological_expression_authority") is not False
    or audit.get("GPL16686_to_GENCODE_v49_crosswalk_complete") is not False
    or audit.get("dependency_install_performed") is not True
    or audit.get("normalization_run_on_biological_arrays") is not False
    or audit.get("ordinary_all_sample_RMA_allowed") is not False
    or audit.get("GEO_series_matrix_allowed_for_model_input") is not False
    or audit.get("model_training_activated") is not False
    or audit.get("sealed_outcomes_read") is not False
):
    raise SystemExit("microarray runtime audit firewall differs")
if not audit["GPL570_expression_summarizer_runtimes"]:
    raise SystemExit("GPL570 expression summarizer inventory is still empty")
if not audit["GPL16686_expression_summarizer_runtimes"]:
    raise SystemExit("GPL16686 expression summarizer inventory is still empty")
freeze_tree(root, {
    "artifact_class": "microarray_activated_runtime_and_inductive_firewall_audit",
    "activation_scope": "single_array_summarizer_and_platform_mapping_runtime_only",
    "R_runtimes_audited": 4,
    "activated_runtime": str(sys.argv[2]),
    "activated_runtime_artifacts_sha256": sys.argv[3],
    "dependency_install_performed": True,
    "capability_probe_arrays_summarized": 2,
    "cohort_expression_matrices_built": False,
    "biological_preprocessing_activated": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed_runtime_activation_with_model_training_still_blocked",
})
verify_frozen_tree(root)
PY
mv "${RUNTIME_STAGE}" "${RUNTIME}"

trap - EXIT
sha256sum "${RUNTIME}/ARTIFACTS.json"
