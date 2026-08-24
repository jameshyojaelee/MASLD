#!/usr/bin/env bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/logs/figures_%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/logs/figures_%j.err

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
PY_BIN="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
cd "$MASLD_PROJECT_ROOT"
"$R_BIN" Analysis/Multimodal_Program_Projection/scripts/05_build_program_context.R
"$R_BIN" scripts/figures/fig5a_input_firewall.R
"$R_BIN" scripts/figures/fig5b_protein_triage.R
"$R_BIN" scripts/figures/composite_mrna_protein.R
"$R_BIN" scripts/figures/fig5d_snatac_accessibility.R
"$PY_BIN" scripts/figures/fig5f_spatial_program_calibration.py
"$R_BIN" scripts/figures/fig5e_multimodal_program_summary.R
"$R_BIN" scripts/figures/update_fig4_manifest.R
"$R_BIN" Analysis/Multimodal_Program_Projection/scripts/06_validate_outputs.R
