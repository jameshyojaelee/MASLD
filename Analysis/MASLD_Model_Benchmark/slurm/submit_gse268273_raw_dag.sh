#!/usr/bin/env bash
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark

if [[ "${GSE268273_SUBMIT:-0}" != 1 ]]; then
  cat <<'COMMANDS'
# Dry-run only. After the declared campaign is authorized, execute with
# GSE268273_SUBMIT=1 bash slurm/submit_gse268273_raw_dag.sh
PLAN_JOB=$(sbatch --parsable slurm/plan_gse268273_raw_campaign_cpu.sbatch)
PLAN_ROOT=".../executions/model-data-075-${PLAN_JOB}/plan"
DOWNLOAD_JOB=$(sbatch --parsable --dependency="afterok:${PLAN_JOB}" --export="ALL,PLAN_ROOT=${PLAN_ROOT}" slurm/download_gse268273_raw_bundles_io.sbatch)
VERIFY_JOB=$(sbatch --parsable --dependency="afterok:${PLAN_JOB}" --export="ALL,PLAN_ROOT=${PLAN_ROOT}" slurm/verify_gse268273_reference_census_cpu.sbatch)
REFERENCE_JOB=$(sbatch --parsable --dependency="afterok:${VERIFY_JOB}" --export="ALL,PLAN_ROOT=${PLAN_ROOT}" slurm/build_gse268273_rsem_reference_cpu.sbatch)
REFERENCE_ROOT=".../executions/model-data-077-${REFERENCE_JOB}/artifact"
QUANT_JOB=$(sbatch --parsable --dependency="afterok:${DOWNLOAD_JOB}:${REFERENCE_JOB}" --export="ALL,PLAN_ROOT=${PLAN_ROOT},REFERENCE_ROOT=${REFERENCE_ROOT}" slurm/quantify_gse268273_rsem_bundles_cpu.sbatch)
FINAL_JOB=$(sbatch --parsable --dependency="afterok:${QUANT_JOB}" --export="ALL,PLAN_ROOT=${PLAN_ROOT},REFERENCE_ROOT=${REFERENCE_ROOT}" slurm/consolidate_gse268273_rsem_counts_cpu.sbatch)
COMMANDS
  exit 0
fi

cd "${ROOT}"
plan_job=$(sbatch --parsable slurm/plan_gse268273_raw_campaign_cpu.sbatch)
plan_job=${plan_job%%;*}
plan_root=${ROOT}/executions/model-data-075-${plan_job}/plan
download_job=$(sbatch \
  --parsable \
  --dependency="afterok:${plan_job}" \
  --export="ALL,PLAN_ROOT=${plan_root}" \
  slurm/download_gse268273_raw_bundles_io.sbatch)
download_job=${download_job%%;*}
verify_job=$(sbatch \
  --parsable \
  --dependency="afterok:${plan_job}" \
  --export="ALL,PLAN_ROOT=${plan_root}" \
  slurm/verify_gse268273_reference_census_cpu.sbatch)
verify_job=${verify_job%%;*}
reference_job=$(sbatch \
  --parsable \
  --dependency="afterok:${verify_job}" \
  --export="ALL,PLAN_ROOT=${plan_root}" \
  slurm/build_gse268273_rsem_reference_cpu.sbatch)
reference_job=${reference_job%%;*}
reference_root=${ROOT}/executions/model-data-077-${reference_job}/artifact
quant_job=$(sbatch \
  --parsable \
  --dependency="afterok:${download_job}:${reference_job}" \
  --export="ALL,PLAN_ROOT=${plan_root},REFERENCE_ROOT=${reference_root}" \
  slurm/quantify_gse268273_rsem_bundles_cpu.sbatch)
quant_job=${quant_job%%;*}
final_job=$(sbatch \
  --parsable \
  --dependency="afterok:${quant_job}" \
  --export="ALL,PLAN_ROOT=${plan_root},REFERENCE_ROOT=${reference_root}" \
  slurm/consolidate_gse268273_rsem_counts_cpu.sbatch)
final_job=${final_job%%;*}
printf 'plan=%s download=%s verify=%s reference=%s quant=%s final=%s\n' \
  "${plan_job}" "${download_job}" "${verify_job}" "${reference_job}" \
  "${quant_job}" "${final_job}"
