#!/bin/bash
#SBATCH --job-name=ag2x_stageB2_analyze
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.err
# Stage B2 analysis. Submit only after the B2 scoring job has ended. No hosted call.
set -euo pipefail
O=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome2/bin/python
export PYTHONHASHSEED=20261111
# New output directory per run; never reuse an existing one.
export AG2X_B2_TAB=$O/stageB2/tables_$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID:-nojob}
if [ -e "$AG2X_B2_TAB" ]; then echo "output dir exists: $AG2X_B2_TAB" >&2; exit 1; fi
mkdir -p $AG2X_B2_TAB
echo "output dir: $AG2X_B2_TAB"
cd $O/code
$PY -m pip freeze > $AG2X_B2_TAB/pip_freeze_analysis.txt
$PY -u b2_20_analyze.py
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $O/code/b2_21_b3_mh.R "$AG2X_B2_TAB"
