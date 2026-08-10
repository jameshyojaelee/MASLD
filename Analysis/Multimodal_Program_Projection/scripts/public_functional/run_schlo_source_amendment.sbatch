#!/bin/bash
#SBATCH --job-name=oshlodata
#SBATCH --partition=io
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/os_h5ad_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/os_h5ad_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJECT_ROOT}"
python3 Analysis/Multimodal_Program_Projection/scripts/public_functional/02c_acquire_schlo_os.py
