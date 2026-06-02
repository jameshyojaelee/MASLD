#!/bin/bash
#SBATCH --job-name=pathway_refactor
#SBATCH --output=logs/pathway_refactor_%j.out
#SBATCH --error=logs/pathway_refactor_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4


# Environment setup previously sourced here - Removed to prevent re-creation
# The environment is already built at the prefix defined below.


echo "Starting Refactored Pipeline Run..."


# Define Execution Command
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
ENV_PREFIX="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/.mamba/pathway_analysis"
RUN_R="$MICROMAMBA run -p $ENV_PREFIX Rscript"

# 1. Run Refactored ssGSEA (Calculates Scores + Validates Genes)
echo "[1/3] Running ssGSEA..."
$RUN_R scripts/run_refactored_ssgsea.R

# 2. Generate Species Contrast Plots
echo "[2/3] Generating Species Contrast Plots..."
$RUN_R scripts/generate_species_contrast_plots.R

# 3. Generate Creative Plots
echo "[3/3] Generating Creative Plots..."
$RUN_R scripts/generate_creative_plots.R


echo "Pipeline Complete!"
