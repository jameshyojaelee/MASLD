#!/usr/bin/env bash
#SBATCH --job-name=p44final
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09/logs/p44final_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09/logs/p44final_%j.err

set -euo pipefail
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CANDIDATE="${ROOT}/Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09"
micromamba run -n rnaseq python \
  "${ROOT}/Analysis/Multimodal_Program_Projection/scripts/multicellular_state_assembly/07_finalize_plan44.py"
micromamba run -n rnaseq Rscript \
  "${ROOT}/Analysis/Multimodal_Program_Projection/scripts/multicellular_state_assembly/08_validate_statistics.R" \
  "${CANDIDATE}"
micromamba run -n rnaseq python \
  "${ROOT}/Analysis/Multimodal_Program_Projection/scripts/multicellular_state_assembly/09_validate_and_freeze.py"
