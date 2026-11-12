#!/bin/bash
#SBATCH --job-name=gset8_fu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/GSET8-20261110T235125Z/logs/slurm-%x-%j.out
set -euo pipefail
bash /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/GSET8-20261110T235125Z/scripts/gset8_followup.sh /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/GSET8-20261110T235125Z/out_23188024 /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/GSET8-20261110T235125Z/out_${SLURM_JOB_ID}
