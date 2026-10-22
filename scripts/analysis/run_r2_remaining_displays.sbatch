#!/bin/bash
#SBATCH --job-name=r2_displays
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=logs/r2_displays_%j.out
#SBATCH --error=logs/r2_displays_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
R2_DISPLAY_ROOT="figures/candidates/r2-current-displays_$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
mkdir -p "$R2_DISPLAY_ROOT/figure5/panels" "$R2_DISPLAY_ROOT/provenance"
echo "$R2_DISPLAY_ROOT"
export MASLD_PROJECT_ROOT="$PWD"
export MPLCONFIGDIR="$R2_DISPLAY_ROOT/provenance/matplotlib"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
R2_PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
R2_RSCRIPT=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
for panel in fig5a_input_firewall fig5b_protein_triage fig5c_mrna_protein_composite fig5d_snatac_accessibility fig5e_multimodal_program_summary fig5f_spatial_program_calibration; do
  cp "figures/main/fig5_molecular_context/panels/$panel.pdf" "$R2_DISPLAY_ROOT/figure5/panels/"
done
sha256sum "$R2_DISPLAY_ROOT"/figure5/panels/*.pdf > "$R2_DISPLAY_ROOT/provenance/figure5_inputs.sha256"
cp scripts/figures/assemble_fig5.py scripts/figures/fig1b_current_interface.R scripts/analysis/synchronize_r2_examples_and_biotypes.py scripts/analysis/run_r2_remaining_displays.sbatch "$R2_DISPLAY_ROOT/provenance/"
"$R2_PYTHON" -m pip freeze > "$R2_DISPLAY_ROOT/provenance/pip_freeze.txt"
"$R2_RSCRIPT" scripts/figures/fig1b_current_interface.R "$R2_DISPLAY_ROOT/figure1"
"$R2_PYTHON" scripts/figures/assemble_fig5.py --candidate-root "$R2_DISPLAY_ROOT"
pdftotext -layout "$R2_DISPLAY_ROOT/figure5/fig5_molecular_context.pdf" "$R2_DISPLAY_ROOT/provenance/figure5_text.txt"
pdftoppm -r 120 -singlefile -png "$R2_DISPLAY_ROOT/figure5/fig5_molecular_context.pdf" "$R2_DISPLAY_ROOT/provenance/figure5_preview"
pdftoppm -r 150 -singlefile -png "$R2_DISPLAY_ROOT/figure1/fig1b_complementary_maps.pdf" "$R2_DISPLAY_ROOT/provenance/figure1b_preview"
"$R2_PYTHON" scripts/analysis/synchronize_r2_examples_and_biotypes.py "$R2_DISPLAY_ROOT/genetics_sources"
pdftoppm -r 150 -singlefile -png "$R2_DISPLAY_ROOT/genetics_sources/Fig2_companion_eqtl_target_gene_biotype.pdf" "$R2_DISPLAY_ROOT/provenance/target_gene_biotype_preview"
sha256sum -c "$R2_DISPLAY_ROOT/provenance/figure5_inputs.sha256" > "$R2_DISPLAY_ROOT/provenance/figure5_copy_check.txt"
printf 'Figure1B and Figure5 vector assembly complete; exact r2 examples and436 source-native target-gene biotypes checked.\n' > "$R2_DISPLAY_ROOT/COMPLETE.txt"
