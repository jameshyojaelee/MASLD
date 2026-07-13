#!/bin/bash
#SBATCH --job-name=metafor
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/heterogeneity_program/phase0/logs/t3_procon_perm_%j.log
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/heterogeneity_program/phase0/logs
export N_PERM=1000 SLURM_CPUS_PER_TASK=16
echo "[start] $(date) host=$(hostname) N_PERM=$N_PERM cores=$SLURM_CPUS_PER_TASK"
~/micromamba/envs/rnaseq/bin/Rscript scripts/heterogeneity_program/phase0/t3_procon.R
echo "[done] $(date)"
