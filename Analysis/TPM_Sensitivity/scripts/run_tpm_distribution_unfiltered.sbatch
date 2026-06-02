#!/usr/bin/env bash
#SBATCH --job-name=tpm_unfiltered
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=cpu
#SBATCH --output=RNA-seq/TPM_analysis/logs/tpm_unfiltered_%j.out
#SBATCH --error=RNA-seq/TPM_analysis/logs/tpm_unfiltered_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p RNA-seq/TPM_analysis/logs

MICROMAMBA=RNA-seq/in-house_MCD_RNAseq/micromamba
ENV_PREFIX=RNA-seq/.mamba/music_deconv
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg

echo "=== TPM Distribution Plots (Unfiltered Raw Counts) ==="
echo "Start time: $(date)"
echo ""

# Debug: Check if file exists
if [ ! -f "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/plot_tpm_distributions_unfiltered.R" ]; then
    echo "Error: R script not found!"
    ls -l /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/
    exit 1
fi

# Set PATH manually to bypass micromamba run issues
export PATH="${ENV_PREFIX}/bin:$PATH"

echo "Rscript location: $(which Rscript)"
echo "R version: $(Rscript --version)"

R --vanilla --no-echo -e "source('/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/plot_tpm_distributions_unfiltered.R'); main()"

echo ""
echo "=== Complete ==="
echo "End time: $(date)"
echo ""
echo "Output files:"
ls -la RNA-seq/plots/tpm_distribution_unfiltered*.pdf 2>/dev/null || echo "No files generated"
