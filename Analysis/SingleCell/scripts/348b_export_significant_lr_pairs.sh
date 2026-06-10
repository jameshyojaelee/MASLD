#!/bin/bash
#SBATCH --job-name=lragg
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/lragg_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/lragg_%j.err

set -eo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "$BASE/Analysis/SingleCell/scripts/logs"
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
Rscript "$BASE/Analysis/SingleCell/scripts/348b_export_significant_lr_pairs.R"
