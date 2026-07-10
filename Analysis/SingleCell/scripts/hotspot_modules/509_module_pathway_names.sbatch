#!/bin/bash
#SBATCH --job-name=msigdbr
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/hotspot_modules/logs/509_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/hotspot_modules/logs/509_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/hotspot_modules/509_module_pathway_names.R
