#!/bin/bash

set -euo pipefail

MODE="${1:-}"
if [[ ! "${MODE}" =~ ^(prepare|run|aggregate)$ ]]; then
  echo "Usage: $0 {prepare|run|aggregate}" >&2
  exit 2
fi

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
FM_ROOT="${PROJECT_ROOT}/GWAS/finemapping"
SCRIPT_ROOT="${FM_ROOT}/src/seqfunc_v2/finemap"
RUN_ROOT="${UNIFORM35_RUN_ROOT:-${FM_ROOT}/runs/uniform35_v2_2026-07-13}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export UNIFORM35_RUN_ROOT="${RUN_ROOT}"
mkdir -p "${RUN_ROOT}/logs"

submit_prepare() {
  local freeze prep finalize
  freeze=$(sbatch --parsable --job-name=SuSiE "${SCRIPT_ROOT}/00_freeze_manifest.sbatch")
  prep=$(sbatch --parsable --job-name=SuSiE --dependency="afterok:${freeze}" \
    --array=1-35%35 "${SCRIPT_ROOT}/01_prepare_study.sbatch")
  finalize=$(sbatch --parsable --job-name=SuSiE --dependency="afterok:${prep}" \
    "${SCRIPT_ROOT}/01b_finalize_preparation.sbatch")
  echo "freeze_job=${freeze}"
  echo "prepare_array_job=${prep}"
  echo "finalize_job=${finalize}"
}

submit_run() {
  local tier task_file n_tasks job partition mem concurrency time_arg
  local -a jobs=()
  [[ -s "${RUN_ROOT}/config/locus_manifest.tsv" ]] || {
    echo "Preparation is not finalized: ${RUN_ROOT}/config/locus_manifest.tsv" >&2
    exit 3
  }
  for tier in small medium large xlarge unknown; do
    task_file="${RUN_ROOT}/config/tasks_${tier}.tsv"
    [[ -f "${task_file}" ]] || continue
    n_tasks=$(( $(wc -l < "${task_file}") - 1 ))
    (( n_tasks > 0 )) || continue
    case "${tier}" in
      # Keep the statistical work single-threaded per locus, but do not impose
      # an artificial array ceiling.  SLURM/QOS availability remains the
      # authoritative concurrency and memory gate.
      small)   partition=cpu;    mem=32G;  concurrency="${n_tasks}"; time_arg=72:00:00 ;;
      medium)  partition=cpu;    mem=64G;  concurrency="${n_tasks}"; time_arg=72:00:00 ;;
      large)   partition=cpu;    mem=128G; concurrency="${n_tasks}"; time_arg=90:00:00 ;;
      xlarge)  partition=bigmem; mem=256G; concurrency="${n_tasks}"; time_arg=90:00:00 ;;
      unknown) partition=bigmem; mem=256G; concurrency="${n_tasks}"; time_arg=90:00:00 ;;
    esac
    job=$(sbatch --parsable --job-name=SuSiE --partition="${partition}" --mem="${mem}" \
      --time="${time_arg}" --array="1-${n_tasks}%${concurrency}" \
      "${SCRIPT_ROOT}/02_run_susie_locus.sbatch" "${task_file}")
    jobs+=("${job}")
    echo "susie_${tier}_job=${job} tasks=${n_tasks} concurrency=${concurrency} mem=${mem} partition=${partition}"
  done
  (( ${#jobs[@]} > 0 )) || { echo "No SuSiE tasks found" >&2; exit 4; }
  local dep aggregate audit
  dep=$(IFS=:; echo "${jobs[*]}")
  aggregate=$(sbatch --parsable --job-name=SuSiE --dependency="afterany:${dep}" "${SCRIPT_ROOT}/03_aggregate.sbatch")
  audit=$(sbatch --parsable --job-name=SuSiE --dependency="afterany:${aggregate}" "${SCRIPT_ROOT}/04_audit.sbatch")
  echo "aggregate_job=${aggregate}"
  echo "audit_job=${audit}"
}

submit_aggregate() {
  local aggregate audit
  aggregate=$(sbatch --parsable --job-name=SuSiE "${SCRIPT_ROOT}/03_aggregate.sbatch")
  audit=$(sbatch --parsable --job-name=SuSiE --dependency="afterany:${aggregate}" "${SCRIPT_ROOT}/04_audit.sbatch")
  echo "aggregate_job=${aggregate}"
  echo "audit_job=${audit}"
}

case "${MODE}" in
  prepare) submit_prepare ;;
  run) submit_run ;;
  aggregate) submit_aggregate ;;
esac
