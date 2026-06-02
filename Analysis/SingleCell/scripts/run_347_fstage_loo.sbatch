#!/usr/bin/env bash
#SBATCH --job-name=347_loo_fstage
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/347_loo_fstage_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/347_loo_fstage_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/347_stage_ccc_loo_dataset.R
