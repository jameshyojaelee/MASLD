#!/bin/bash
#SBATCH --job-name=eqtl_psd_E
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/logs/overlap_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/logs/overlap_%j.err
set -eo pipefail
export PYTHONNOUSERSITE=1
export RUN_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY ${RUN_DIR}/code/overlap_refit.py ${RUN_DIR}/results/overlap_$SLURM_JOB_ID
$PY -m pip freeze > ${RUN_DIR}/results/overlap_$SLURM_JOB_ID/pip_freeze.txt
