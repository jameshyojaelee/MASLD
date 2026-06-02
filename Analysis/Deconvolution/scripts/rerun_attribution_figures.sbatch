#!/usr/bin/env bash
#SBATCH --job-name=attrib_figures
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=Analysis/Deconvolution/logs/attrib_figures_%j.out
#SBATCH --error=Analysis/Deconvolution/logs/attrib_figures_%j.err

# Re-run deconvolution attribution (script 25) + BayesPrism validation (script 22)
# + regenerate figure 3 (deconvolution only; single-cell panels moved to fig5)

set -euo pipefail

echo "=== Attribution + Figure Re-generation ==="
echo "Date: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export R_LIBS_USER=""

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Step 1: Re-run deconvolution attribution (script 25)
echo ""
echo "=== Script 25: Deconvolution Attribution ==="
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/25_deconv_attribution.R

# Step 2: Re-run BayesPrism validation (script 22)
echo ""
echo "=== Script 22: BayesPrism Cross-Validation ==="
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/22_bayesprism_deconvolution.R || \
    echo "WARNING: Script 22 failed (non-fatal)"

# Step 3: Regenerate figures
echo ""
echo "=== Figure 3: Deconvolution Hero ==="
Rscript scripts/figures/fig3_deconvolution_hero.R

echo ""
echo "=== Done ==="
echo "End: $(date)"
ls -la figures/fig3_deconvolution_hero.pdf 2>/dev/null
