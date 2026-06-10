#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/_c2_rerender/logs/sweep_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/_c2_rerender/logs/sweep_%j.err
# Driven by run_one.sh over a chunk file passed as $1.
set -uo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
bash "$BASE/scripts/figures/_c2_rerender/run_one.sh" "$1"
