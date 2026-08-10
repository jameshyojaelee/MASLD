#!/usr/bin/env bash
#SBATCH --job-name=p44human
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09/logs/p44human_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09/logs/p44human_%j.err

set -euo pipefail
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CANDIDATE="${ROOT}/Analysis/Multimodal_Program_Projection/candidates/multicellular-assembly-response-2026-08-09"
micromamba run -n rnaseq Rscript \
  "${ROOT}/Analysis/Multimodal_Program_Projection/scripts/multicellular_state_assembly/06_analyze_gse175448_human.R" \
  "${CANDIDATE}"
