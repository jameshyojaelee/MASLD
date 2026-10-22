#!/bin/bash
#SBATCH --job-name=N3amb_decontx
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/N3-ambient-20261013T060951Z/%x-%j.out
# N3 ambient check, job A: per-library decontX with empty-droplet profiles.
# Usage: sbatch run_decontx.sbatch [MAX_LIBRARIES]   (no argument = all libraries)
set -euo pipefail
export REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export RES=$REPO/Analysis/Hypotheses/N3-ambient-20261013T060951Z
export SCR=$RES/scripts
N3=$REPO/Analysis/Hypotheses/N3-20261009T140326Z
export PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export RS_DECONTX=/gpfs/commons/home/jameslee/micromamba/envs/decontx/bin/Rscript
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
MAXLIB="${1:-0}"
test -f "$RES/PRESPEC.md"
# node-local temporary files of this job only (one library per worker at a time)
export WORK=$(mktemp -d "${TMPDIR:-/tmp}/n3amb_${SLURM_JOB_ID}_XXXX")
trap 'rm -rf "$WORK"' EXIT
echo "host=$(hostname) job=$SLURM_JOB_ID work=$WORK start=$(date -u)"
(cd "$SCR" && sha256sum *.py *.R *.sh *.sbatch) > "$RES/script_sha256_job${SLURM_JOB_ID}.txt"
sha256sum "$RES/PRESPEC.md" >> "$RES/script_sha256_job${SLURM_JOB_ID}.txt"
"$PY" -I "$SCR/00_lists.py" "$REPO" "$N3" "$RES"
LIST="$RES/library_list.tsv"
if [ "$MAXLIB" -gt 0 ]; then head -n "$MAXLIB" "$LIST" > "$WORK/list.tsv"; else cp "$LIST" "$WORK/list.tsv"; fi
xargs -P 4 -L 1 bash "$SCR/run_one_library.sh" < "$WORK/list.tsv"
echo "per-library results: $(ls "$RES/per_library" | wc -l) end=$(date -u)"
