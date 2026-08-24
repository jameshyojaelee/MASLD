#!/usr/bin/env bash
#SBATCH --job-name=fig4e_promote
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4e_promote_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4e_promote_%j.err

set -euo pipefail
python3 /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/validate_fig4e_uniform_promotion.py
