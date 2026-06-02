#!/usr/bin/env bash
#SBATCH --job-name=348_bulkconc_fstage
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/348_bulkconc_fstage_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/348_bulkconc_fstage_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/348_stage_ccc_bulk_concordance.R
