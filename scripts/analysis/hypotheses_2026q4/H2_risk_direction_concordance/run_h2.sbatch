#!/bin/bash
#SBATCH --job-name=H2_riskdir
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
# H2 risk-direction concordance. Usage: sbatch --output=<result_dir>/logs/%j_h2.out run_h2.sbatch <result_dir>
set -euo pipefail
out=$1
here=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H2_risk_direction_concordance
py=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
cd /
echo "host $(hostname) job ${SLURM_JOB_ID:-none} start $(date -u +%FT%TZ)"
test -f "$out/PRESPEC_H2.md"
run=$out/run_${SLURM_JOB_ID}
mkdir "$run"
echo 20260923 > "$run/SEED.txt"
$py -I -m pip freeze > "$run/pip_freeze.txt"
sha256sum "$here"/*.py "$here"/*.sbatch "$out/PRESPEC_H2.md" > "$run/code_sha256.txt"
$py -I "$here/h2_risk_direction.py" "$out" "$run"
echo "end $(date -u +%FT%TZ)"
