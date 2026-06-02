#!/bin/bash
#SBATCH --job-name=b6_breakpoint
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/primary/scripts/b6_breakpoint_%j.log

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/primary/scripts
Rscript b6_breakpoint_sensitivity.R
