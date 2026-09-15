#!/usr/bin/env bash
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export PYTHONPATH="${ROOT}/src:${ROOT}"

LOCK_FILES=(
  scripts/plan_gse105127_production_dag.py
  scripts/freeze_gse105127_reference_bundle.py
  scripts/download_gse105127_rna_bundle.py
  scripts/collapse_gse105127_rrbs_bundle.py
  scripts/build_gse105127_rsem_reference.py
  scripts/quantify_gse105127_rna_bundle.py
  scripts/consolidate_gse105127_rna_quantification.py
  scripts/build_gse105127_cpg_crosswalk.py
  scripts/finalize_gse105127_production_activation.py
  tests/unit/test_build_gse105127_assay_native_activation.py
  tests/unit/test_gse105127_production_dag.py
  config/evaluation/gse105127_adjacent_section_rna_rrbs_task.toml
  slurm/gse105127_production_common.sh
  slurm/gse105127_stage_plan_cpu.sbatch
  slurm/gse105127_stage_reference_io.sbatch
  slurm/gse105127_stage_rrbs_io.sbatch
  slurm/gse105127_stage_rna_download_io.sbatch
  slurm/gse105127_stage_rsem_reference_cpu.sbatch
  slurm/gse105127_stage_crosswalk_io.sbatch
  slurm/gse105127_stage_rna_quant_cpu.sbatch
  slurm/gse105127_stage_consolidate_cpu.sbatch
  slurm/gse105127_stage_finalize_cpu.sbatch
  slurm/validate_gse105127_production_dag_cpu.sbatch
  slurm/submit_gse105127_production_dag.sh
  executions/gse105127-assay-native-activation-21082444/ARTIFACTS.json
  executions/gse105127-rrbs-source-qc-r2-21066008/ARTIFACTS.json
  executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json
)
LOCK_PAYLOAD=$(sha256sum "${LOCK_FILES[@]}")
CAMPAIGN_HASH=$(printf '%s\n' "${LOCK_PAYLOAD}" | sha256sum | cut -d' ' -f1)
GSE105127_CAMPAIGN_ID=${CAMPAIGN_HASH:0:16}
export GSE105127_CAMPAIGN_ID
CAMPAIGN=${ROOT}/executions/gse105127-production-${GSE105127_CAMPAIGN_ID}
CONTROL=${CAMPAIGN}/control
MINIMUM_FREE_BYTES=$((500 * 1024 * 1024 * 1024))
AVAILABLE_BYTES=$(df -PB1 "${ROOT}/executions" | awk 'NR == 2 {print $4}')
if [[ ! ${AVAILABLE_BYTES} =~ ^[0-9]+$ || ${AVAILABLE_BYTES} -lt ${MINIMUM_FREE_BYTES} ]]; then
  printf 'Production DAG requires at least 500 GiB filesystem free space\n' >&2
  exit 4
fi

if [[ -f ${CONTROL}/COMPLETE ]]; then
  "${PYTHON}" - "${CONTROL}" "${CAMPAIGN_HASH}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1])
verify_frozen_tree(root)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
if receipt.get("campaign_hash") != sys.argv[2]:
    raise SystemExit("campaign control hash differs")
PY
else
  if [[ -e ${CONTROL} ]]; then
    printf 'Incomplete campaign control requires manual audit: %s\n' "${CONTROL}" >&2
    exit 2
  fi
  mkdir -p "${CAMPAIGN}"
  CONTROL_STAGE=${CAMPAIGN}/.control.staging
  if [[ -e ${CONTROL_STAGE} ]]; then
    printf 'Refusing to replace campaign control stage\n' >&2
    exit 3
  fi
  mkdir "${CONTROL_STAGE}"
  printf '%s\n' "${LOCK_PAYLOAD}" > "${CONTROL_STAGE}/source-lock.sha256"
  "${PYTHON}" - "${CONTROL_STAGE}" "${CONTROL}" "${CAMPAIGN_HASH}" "${GSE105127_CAMPAIGN_ID}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, publish_directory_noreplace, verify_frozen_tree

stage, final = map(Path, sys.argv[1:3])
receipt = {
    "schema_version": "masld-bench-gse105127-production-control-v1",
    "status": "source_locked_before_submission",
    "campaign_hash": sys.argv[3],
    "campaign_id": sys.argv[4],
    "participants": 19,
    "participant_zone_rows": 57,
    "bundles": 8,
    "production_slurm_tasks": 30,
    "production_cpu_hour_ceiling": 5336,
    "production_memory_gb_hour_ceiling": 44192,
    "gpu_hours": 0,
    "persistent_storage_ceiling_gib": 200,
    "transient_storage_ceiling_gib": 300,
    "labels_accessed": False,
    "fit_or_score_performed": False,
}
(stage / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
freeze_tree(stage, {"artifact_class": "gse105127_production_control", "status": "locked"})
verify_frozen_tree(stage)
publish_directory_noreplace(stage, final)
verify_frozen_tree(final)
PY
fi
mkdir -p "${CAMPAIGN}/attempts" "${CAMPAIGN}/submissions" executions/logs
(
  cd "${ROOT}"
  sha256sum --check --strict "${CONTROL}/source-lock.sha256"
)

submit_job() {
  local job_id
  job_id=$(sbatch --parsable --export=ALL,GSE105127_CAMPAIGN_ID="${GSE105127_CAMPAIGN_ID}" "$@")
  printf '%s\n' "${job_id%%;*}"
}

VALIDATE=$(submit_job slurm/validate_gse105127_production_dag_cpu.sbatch)
PLAN=$(submit_job --dependency="afterok:${VALIDATE}" slurm/gse105127_stage_plan_cpu.sbatch)
REFERENCE=$(submit_job --dependency="afterok:${PLAN}" slurm/gse105127_stage_reference_io.sbatch)
RRBS=$(submit_job --dependency="afterok:${PLAN}" slurm/gse105127_stage_rrbs_io.sbatch)
RNA_DOWNLOAD=$(submit_job --dependency="afterok:${PLAN}" slurm/gse105127_stage_rna_download_io.sbatch)
RSEM_REFERENCE=$(submit_job --dependency="afterok:${PLAN}" slurm/gse105127_stage_rsem_reference_cpu.sbatch)
CROSSWALK=$(submit_job --dependency="afterok:${REFERENCE}:${RRBS}" slurm/gse105127_stage_crosswalk_io.sbatch)
RNA_QUANT=$(submit_job --dependency="afterok:${RNA_DOWNLOAD}:${RSEM_REFERENCE}" slurm/gse105127_stage_rna_quant_cpu.sbatch)
CONSOLIDATE=$(submit_job --dependency="afterok:${RNA_QUANT}" slurm/gse105127_stage_consolidate_cpu.sbatch)
FINALIZE=$(submit_job --dependency="afterok:${CROSSWALK}:${CONSOLIDATE}" slurm/gse105127_stage_finalize_cpu.sbatch)

SUBMISSION_STAGE=${CAMPAIGN}/submissions/.submission-${VALIDATE}.staging
SUBMISSION=${CAMPAIGN}/submissions/submission-${VALIDATE}
mkdir "${SUBMISSION_STAGE}"
printf 'stage\tjob_id\tdependency\nvalidation\t%s\t\nplan\t%s\tafterok:%s\nreference\t%s\tafterok:%s\nrrbs_array\t%s\tafterok:%s\nrna_download_array\t%s\tafterok:%s\nrsem_reference\t%s\tafterok:%s\ncrosswalk\t%s\tafterok:%s:%s\nrna_quant_array\t%s\tafterok:%s:%s\nconsolidate\t%s\tafterok:%s\nfinalize\t%s\tafterok:%s:%s\n' \
  "${VALIDATE}" "${PLAN}" "${VALIDATE}" "${REFERENCE}" "${PLAN}" \
  "${RRBS}" "${PLAN}" "${RNA_DOWNLOAD}" "${PLAN}" "${RSEM_REFERENCE}" "${PLAN}" \
  "${CROSSWALK}" "${REFERENCE}" "${RRBS}" "${RNA_QUANT}" "${RNA_DOWNLOAD}" "${RSEM_REFERENCE}" \
  "${CONSOLIDATE}" "${RNA_QUANT}" "${FINALIZE}" "${CROSSWALK}" "${CONSOLIDATE}" \
  > "${SUBMISSION_STAGE}/jobs.tsv"
"${PYTHON}" - "${SUBMISSION_STAGE}" "${SUBMISSION}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, publish_directory_noreplace, verify_frozen_tree

stage, final = map(Path, sys.argv[1:])
freeze_tree(stage, {"artifact_class": "gse105127_production_submission", "status": "submitted"})
verify_frozen_tree(stage)
publish_directory_noreplace(stage, final)
verify_frozen_tree(final)
PY
printf 'campaign=%s\nsubmission=%s\nfinal_job=%s\n' "${CAMPAIGN}" "${SUBMISSION}" "${FINALIZE}"
