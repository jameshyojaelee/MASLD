#!/bin/bash
#SBATCH --job-name=teamF2_figS_regen
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --output=scripts/figures/logs/teamF2_figS_regen_%j.out
#SBATCH --error=scripts/figures/logs/teamF2_figS_regen_%j.err
set -o pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

mkdir -p scripts/figures/logs

for s in figS_sex_dimorphism.R figS_network_v4.R figS_healthy_control_audit.R figS_transition_tf_activity.R; do
  echo "[$(date)] running $s"
  Rscript scripts/figures/$s 2>&1 | tee scripts/figures/logs/${s%.R}_regen.log
done

echo "[$(date)] team F2 regen complete"
