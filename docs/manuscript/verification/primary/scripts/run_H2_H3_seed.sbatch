#!/bin/bash
#SBATCH --job-name=H2H3_seed
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/primary/scripts/H2H3_seed_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/primary/scripts/H2H3_seed_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/verification/primary/scripts

# Start with 200 draws (fast); update with 1000 when scaling up
N_DRAWS=${N_DRAWS:-200} python H2_H3_independent_seed.py

echo "Job completed at $(date)"
