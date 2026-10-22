#!/bin/bash
#SBATCH --job-name=H17_terminal_check
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=logs/H17_terminal_check_%j.out
#SBATCH --error=logs/H17_terminal_check_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export MASLD_PROJECT_ROOT="$PWD"
export H17_CHECK_ROOT="$PWD/Analysis/Hypotheses/H17-terminal-check-20261014_${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python scripts/analysis/hypotheses_2026q4/H17_lipid_response_core/check_h17_source_corrected.py
cp scripts/analysis/hypotheses_2026q4/H17_lipid_response_core/check_h17_source_corrected.py scripts/analysis/hypotheses_2026q4/H17_lipid_response_core/run_check_h17_source_corrected.sbatch "$H17_CHECK_ROOT/"
