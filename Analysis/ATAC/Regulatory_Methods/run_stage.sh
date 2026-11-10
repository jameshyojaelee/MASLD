#!/usr/bin/env bash
# Run one candidate stage inside a verified compute allocation. Never submits jobs.
set -euo pipefail
TASK_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CONFIG=${REGULATORY_CONFIG:-"${TASK_ROOT}/config/analysis.json"}
SNAP_PY=/gpfs/commons/home/jameslee/micromamba/envs/snapatac2/bin/python
MODEL_PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
RNA_R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ -z ${SLURM_JOB_ID:-} ]]; then
  echo 'Run analysis stages inside a SLURM compute allocation.' >&2
  exit 2
fi
NODE_RECORD=$(scontrol show job "$SLURM_JOB_ID" -o)
ALLOCATED_NODES=$(sed -n 's/.* NodeList=\([^ ]*\).*/\1/p' <<< "$NODE_RECORD")
if ! scontrol show hostnames "$ALLOCATED_NODES" | grep -Fxq "$(hostname -s)"; then
  echo 'Current host is not in the recorded SLURM allocation.' >&2
  exit 2
fi
export SLURM_JOB_NODELIST="$ALLOCATED_NODES"
STAGE=${1:?Usage: run_stage.sh STAGE OUTPUT [stage arguments]}
OUTPUT=${2:?A new output location is required}
shift 2
case "$STAGE" in
  audit)
    exec "$SNAP_PY" "$TASK_ROOT/scripts/input_counts.py" audit --config "$CONFIG" --output-root "$OUTPUT" "$@"
    ;;
  count-pilot)
    exec "$SNAP_PY" "$TASK_ROOT/scripts/input_counts.py" pilot --config "$CONFIG" --output-root "$OUTPUT" "$@"
    ;;
  count)
    for COHORT in GSE244832 GSE281367; do
      for LINEAGE in hepatocyte stellate macrophage; do
        "$SNAP_PY" "$TASK_ROOT/scripts/input_counts.py" count --config "$CONFIG" --output-root "$OUTPUT" --cohort "$COHORT" --lineage "$LINEAGE" "$@"
      done
    done
    ;;
  da)
    exec "$RNA_R" "$TASK_ROOT/scripts/compare_da.R" --config "$CONFIG" --output "$OUTPUT" "$@"
    ;;
  paired-pilot|paired)
    MODE=full
    [[ "$STAGE" == paired-pilot ]] && MODE=pilot
    "$RNA_R" "$TASK_ROOT/scripts/export_paired.R" --config "$CONFIG" --output "$OUTPUT" --mode "$MODE" --lineage all "$@"
    exec "$SNAP_PY" "$TASK_ROOT/scripts/qualify_sources.py" export-rna --config "$CONFIG" --paired-root "$OUTPUT"
    ;;
  network)
    exec "$MODEL_PY" "$TASK_ROOT/scripts/network_compare.py" --config "$CONFIG" --output "$OUTPUT" "$@"
    ;;
  network-sensitivity)
    exec "$MODEL_PY" "$TASK_ROOT/scripts/network_sensitivity.py" --config "$CONFIG" --output "$OUTPUT" "$@"
    ;;
  reference)
    exec "$MODEL_PY" "$TASK_ROOT/scripts/evaluate_link_reference.py" --output "$OUTPUT" "$@"
    ;;
  motif-context)
    exec "$MODEL_PY" "$TASK_ROOT/scripts/annotate_motif_context.py" --output "$OUTPUT" "$@"
    ;;
  summarize)
    exec "$MODEL_PY" "$TASK_ROOT/scripts/summarize_candidate.py" --output "$OUTPUT" "$@"
    ;;
  figures)
    exec "$RNA_R" "$TASK_ROOT/scripts/render_results.R" --output "$OUTPUT" "$@"
    ;;
  *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
