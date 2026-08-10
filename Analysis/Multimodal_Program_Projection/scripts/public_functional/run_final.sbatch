#!/bin/bash
#SBATCH --job-name=pubfinal
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/final_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/final_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJECT_ROOT}"
micromamba run -n rapids_singlecell python Analysis/Multimodal_Program_Projection/scripts/public_functional/16_capture_environment.py
micromamba run -n rapids_singlecell python Analysis/Multimodal_Program_Projection/scripts/public_functional/11_run_null_permutations.py
micromamba run -n rapids_singlecell python Analysis/Multimodal_Program_Projection/scripts/public_functional/12_build_final_integration.py
micromamba run -n rapids_singlecell python Analysis/Multimodal_Program_Projection/scripts/public_functional/13_validate_final.py
