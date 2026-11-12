#!/bin/bash
#SBATCH --job-name=ag2x_stageB_b3a_fix
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.err
set -euo pipefail
O=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome2/bin/python
export PYTHONHASHSEED=20261111
mkdir -p $O/stageB/tables
cd $O/code
$PY -u b23_b3a_corrected.py
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $O/code/b24_b3a_corrected_mh.R
