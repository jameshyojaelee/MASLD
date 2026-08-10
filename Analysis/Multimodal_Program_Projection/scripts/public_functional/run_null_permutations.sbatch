#!/bin/bash
#SBATCH --job-name=pubnulls
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/nulls_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/nulls_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJECT_ROOT}"
micromamba run -n rapids_singlecell python Analysis/Multimodal_Program_Projection/scripts/public_functional/11_run_null_permutations.py
