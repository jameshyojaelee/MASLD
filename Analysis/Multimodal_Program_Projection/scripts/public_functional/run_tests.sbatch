#!/bin/bash
#SBATCH --job-name=pubtests
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/tests_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/logs/tests_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJECT_ROOT}"
python3 -m py_compile Analysis/Multimodal_Program_Projection/scripts/public_functional/*.py
micromamba run -n rnaseq Rscript Analysis/Multimodal_Program_Projection/scripts/public_functional/test_public_functional.R
