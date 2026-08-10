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
MAX_ARRAY_CONCURRENCY="${UNIFORM35_MAX_CONCURRENCY:-200}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export UNIFORM35_RUN_ROOT="${RUN_ROOT}"
mkdir -p "${RUN_ROOT}/logs"

# The sbatch files declare RELATIVE log paths, resolved against the submission
# working directory; cd here so they land in this run's logs/ rather than
# whichever run was hardcoded.  Every submission below ALSO passes --output and
# --error explicitly, which overrides the directive outright, and sub() then
# asserts SLURM actually resolved them inside RUN_ROOT.
cd "${RUN_ROOT}"
echo "run_root=${RUN_ROOT}"

# sbatch wrapper: force this run's log paths and verify SLURM honoured them.
# A job whose StdOut escapes RUN_ROOT is cancelled rather than left to
# silently overwrite another run's logs.
sub() {
  local stem="$1"; shift
  local jid
  jid=$(sbatch --parsable --job-name=SuSiE \
    --output="${RUN_ROOT}/logs/${stem}.out" \
    --error="${RUN_ROOT}/logs/${stem}.err" "$@")
  local stdout_path
  stdout_path=$(scontrol show job "${jid}" 2>/dev/null | tr ' ' '\n' | grep '^StdOut=' | head -1 | cut -d= -f2-)
  if [[ -n "${stdout_path}" && "${stdout_path}" != "${RUN_ROOT}"* ]]; then
    echo "FATAL: job ${jid} would write logs outside the run root: ${stdout_path}" >&2
    scancel "${jid}" 2>/dev/null || true
    exit 5
  fi
  echo "${jid}"
}

submit_prepare() {
  local freeze prep finalize
  freeze=$(sub "freeze_%j" "${SCRIPT_ROOT}/00_freeze_manifest.sbatch")
  prep=$(sub "prep_%A_%a" --dependency="afterok:${freeze}" \
    --array=1-35%35 "${SCRIPT_ROOT}/01_prepare_study.sbatch")
  finalize=$(sub "finalize_%j" --dependency="afterok:${prep}" \
    "${SCRIPT_ROOT}/01b_finalize_preparation.sbatch")
  echo "freeze_job=${freeze}"
  echo "prepare_array_job=${prep}"
  echo "finalize_job=${finalize}"
}

submit_run() {
  local tier task_file n_tasks job partition mem cpus concurrency time_arg
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
    # Memory is requested against the modelled peak (~48*M^2 bytes; see
    # predicted_peak_bytes in common.R), not the old flat ladder.  The largest
    # block in the portfolio (M = 23,331) predicts ~26 GB, so nothing needs
    # bigmem -- the previous xlarge/unknown tiers requested 256G on bigmem for
    # loci with a real footprint around 14 GB and queued badly for no benefit.
    # large/xlarge get 4 CPUs because the ridge ladder runs up to 7 Choleskys,
    # which are BLAS-3 and thread well.
    case "${tier}" in
      small)   partition=cpu; mem=12G;  cpus=1; time_arg=72:00:00 ;;
      medium)  partition=cpu; mem=32G;  cpus=1; time_arg=72:00:00 ;;
      large)   partition=cpu; mem=64G;  cpus=4; time_arg=90:00:00 ;;
      xlarge)  partition=cpu; mem=128G; cpus=4; time_arg=90:00:00 ;;
      unknown) partition=cpu; mem=64G;  cpus=4; time_arg=90:00:00 ;;
    esac
    # Cap concurrency: ~2.2 TB of LD reads at 1,950-wide would thrash GPFS for
    # every other user of the filesystem.
    concurrency=$(( n_tasks < MAX_ARRAY_CONCURRENCY ? n_tasks : MAX_ARRAY_CONCURRENCY ))
    job=$(sub "susie_%A_%a" --partition="${partition}" --mem="${mem}" \
      --cpus-per-task="${cpus}" \
      --time="${time_arg}" --array="1-${n_tasks}%${concurrency}" \
      "${SCRIPT_ROOT}/02_run_susie_locus.sbatch" "${task_file}")
    jobs+=("${job}")
    echo "susie_${tier}_job=${job} tasks=${n_tasks} concurrency=${concurrency} mem=${mem} cpus=${cpus} partition=${partition}"
  done
  (( ${#jobs[@]} > 0 )) || { echo "No SuSiE tasks found" >&2; exit 4; }
  local dep aggregate audit
  dep=$(IFS=:; echo "${jobs[*]}")
  aggregate=$(sub "aggregate_%j" --dependency="afterany:${dep}" "${SCRIPT_ROOT}/03_aggregate.sbatch")
  audit=$(sub "audit_%j" --dependency="afterany:${aggregate}" "${SCRIPT_ROOT}/04_audit.sbatch")
  echo "aggregate_job=${aggregate}"
  echo "audit_job=${audit}"
}

submit_aggregate() {
  local aggregate audit
  aggregate=$(sub "aggregate_%j" "${SCRIPT_ROOT}/03_aggregate.sbatch")
  audit=$(sub "audit_%j" --dependency="afterany:${aggregate}" "${SCRIPT_ROOT}/04_audit.sbatch")
  echo "aggregate_job=${aggregate}"
  echo "audit_job=${audit}"
}

case "${MODE}" in
  prepare) submit_prepare ;;
  run) submit_run ;;
  aggregate) submit_aggregate ;;
esac
