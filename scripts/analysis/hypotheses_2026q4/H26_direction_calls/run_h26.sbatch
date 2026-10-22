#!/bin/bash
#SBATCH --job-name=H26_direction
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=90:00:00
# H26 direction calls. Usage: sbatch --output=<result_dir>/logs/%j_h26.out run_h26.sbatch <result_dir>
set -euo pipefail
out=$1
here=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H26_direction_calls
py=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
cd /
echo "host $(hostname) job ${SLURM_JOB_ID:-none} start $(date -u +%FT%TZ)"
test -f "$out/PRESPEC_H26.md"
run=$out/run_${SLURM_JOB_ID}
mkdir "$run"
$py -I -m pip freeze > "$run/pip_freeze.txt"
sha256sum "$here"/*.py "$here"/*.sbatch "$out/PRESPEC_H26.md" > "$run/code_sha256.txt"
$py -I "$here/h26_direction_calls.py" "$out" "$run"
echo "end $(date -u +%FT%TZ)"
