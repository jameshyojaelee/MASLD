#!/bin/bash
#SBATCH --job-name=masld_review_read_phase
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --dependency=afterok:21773997
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_read_phase_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_read_phase_%j.err
set -euo pipefail
umask 077
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
review_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
review_env=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
review_base=$review_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
"$review_env/python" "$review_root/scripts/analysis/alphagenome_campaign/review_phase_counts.py" \
 --out "$review_base/review/phase-counts-$SLURM_JOB_ID"
"$review_env/Rscript" "$review_root/scripts/analysis/alphagenome_campaign/review_currin_crosswalk.R" \
 "$review_root" "$review_base/review/currin-crosswalk-$SLURM_JOB_ID"
