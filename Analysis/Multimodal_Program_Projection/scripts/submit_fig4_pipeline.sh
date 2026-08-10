#!/usr/bin/env bash

set -euo pipefail
ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$ROOT"
mkdir -p Analysis/Multimodal_Program_Projection/logs

SBATCH_OVERRIDE=()
if [[ -n "${FIG4_PARTITION:-}" ]]; then SBATCH_OVERRIDE+=("--partition=${FIG4_PARTITION}"); fi
if [[ -n "${FIG4_QOS:-}" ]]; then SBATCH_OVERRIDE+=("--qos=${FIG4_QOS}"); fi

PROTEIN_JOB=$(sbatch --parsable "${SBATCH_OVERRIDE[@]}" Analysis/Multimodal_Program_Projection/scripts/run_registry_proteomics.sbatch)
ATAC_JOB=$(sbatch --parsable "${SBATCH_OVERRIDE[@]}" --dependency="afterok:${PROTEIN_JOB}" Analysis/Multimodal_Program_Projection/scripts/run_atac.sbatch)
SPATIAL_JOB=$(sbatch --parsable "${SBATCH_OVERRIDE[@]}" --dependency="afterok:${PROTEIN_JOB}" Analysis/Multimodal_Program_Projection/scripts/run_spatial.sbatch)
FIGURE_JOB=$(sbatch --parsable "${SBATCH_OVERRIDE[@]}" --dependency="afterok:${ATAC_JOB}:${SPATIAL_JOB}" Analysis/Multimodal_Program_Projection/scripts/run_figures.sbatch)

printf 'proteomics\t%s\nATAC\t%s\nspatial\t%s\nfigures\t%s\n' \
  "$PROTEIN_JOB" "$ATAC_JOB" "$SPATIAL_JOB" "$FIGURE_JOB"
