#!/bin/bash
#SBATCH --job-name=H2H3_quick
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=8:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/primary/scripts/H2H3_quick_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/primary/scripts/H2H3_quick_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/primary/scripts

N_DRAWS=50 python H2_H3_independent_seed.py
