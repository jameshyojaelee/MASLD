#!/bin/bash
#SBATCH --job-name=ag2x_stageC_score
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.err
# Stage C scoring: sends anonymous count matrices to the hosted API. Submit only after the orchestrator's go-ahead.
set -euo pipefail
O=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome2/bin/python
export PYTHONHASHSEED=20261111 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
mkdir -p $O/stageC/raw $O/stageC/tables
cd $O/code
$PY -m pip freeze > $O/stageC/raw/pip_freeze_scoring.txt
sha256sum $O/PRESPEC_AG2X_stageC.md $O/code/c01_stageC_lib.py $O/code/c10_score_stageC.py $O/code/c20_analyze_stageC.py > $O/stageC/raw/c10_script_sha256_${SLURM_JOB_ID}.txt
$PY -u c10_score_stageC.py "$@"
if [ ! -f $O/stageC/NOT_TESTABLE.json ] && [ -f $O/stageC/raw/scoring_summary.json ]; then
  $PY -u c20_analyze_stageC.py
  $PY -u c20_analyze_stageC.py --all-tracks
fi
