#!/usr/bin/env bash
#SBATCH --job-name=model-data-073
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PREPARATION=${ROOT}/executions/model-data-072-21080651-preparation
PREPARATION_ARTIFACTS_SHA256=8dd6f7288f7c06d7e1db5e64d87eedb3cfa122c13afee3da6293ff204144f96b
INDUCTIVE=${ROOT}/executions/model-data-072-21080651-inductive
INDUCTIVE_ARTIFACTS_SHA256=ea27d2ae95b914b1e01831798a37db692f681d504f6583606f0c902439a0cd58
PREPROCESSING=${ROOT}/config/evaluation/microarray_transfer_preprocessing.toml
PREPROCESSING_SHA256=11162c381017bd16736433f478c0c82d97726ae2d896957d27efa5aee899c27e
LABELS=${ROOT}/config/evaluation/microarray_transfer_label_semantics.toml
LABELS_SHA256=b66c20ef5c6bc2ce7597affe7816358bd7fd4ca65df8e3f4be1096dbddaf0513

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
  | sha256sum --check --strict

module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python - "${PREPARATION}" "${INDUCTIVE}" <<'PY'
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

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python scripts/finalize_microarray_runtime_audit.py \
  --preparation "${PREPARATION}" \
  --inductive-fixture "${INDUCTIVE}" \
  --runtime-receipt "${RECEIPTS_STAGE}/module_R_4.4.3.tsv" \
  --runtime-receipt "${RECEIPTS_STAGE}/module_R_4.4.1.tsv" \
  --runtime-receipt "${RECEIPTS_STAGE}/micromamba_rnaseq.tsv" \
  --preprocessing-contract "${PREPROCESSING}" \
  --label-contract "${LABELS}" \
  --output "${RUNTIME_STAGE}"
mv "${RECEIPTS_STAGE}"/*.tsv "${RUNTIME_STAGE}/"
rmdir "${RECEIPTS_STAGE}"
sha256sum \
  scripts/audit_microarray_r_runtime.R \
  scripts/finalize_microarray_runtime_audit.py \
  tests/unit/test_finalize_microarray_runtime_audit.py \
  "${PREPROCESSING}" "${LABELS}" \
  "${PREPARATION}/ARTIFACTS.json" "${INDUCTIVE}/ARTIFACTS.json" \
  > "${RUNTIME_STAGE}/source.sha256"
python - "${RUNTIME_STAGE}" <<'PY'
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
    or audit.get("dependency_install_performed") is not False
    or audit.get("model_training_activated") is not False
):
    raise SystemExit("microarray runtime audit firewall differs")
freeze_tree(root, {
    "artifact_class": "microarray_existing_runtime_and_inductive_firewall_audit",
    "retry_scope": "runtime_receipt_only_reusing_verified_frozen_inputs",
    "R_runtimes_audited": 3,
    "dependency_install_performed": False,
    "biological_preprocessing_activated": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed_runtime_inventory_biological_activation_blocked",
})
verify_frozen_tree(root)
PY
mv "${RUNTIME_STAGE}" "${RUNTIME}"

trap - EXIT
sha256sum "${RUNTIME}/ARTIFACTS.json"
