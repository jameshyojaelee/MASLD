#!/bin/bash
#SBATCH --job-name=lfc_sweep
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/lfc_sweep_%j.out
#SBATCH --error=logs/lfc_sweep_%j.err

# B3 LFC × padj grid sweep — reuses pre-computed LOO dream results.
# I/O light: reads 6 CSVs (~35K rows each), no dream re-run.

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/RNA-seq/scripts/lfc_sweep

source /gpfs/commons/home/jameslee/.bashrc
micromamba activate rnaseq
set -u

Rscript lfc_grid_sweep.R
echo "Sweep done. Plotting..."
Rscript plot_lfc_sweep.R
echo "All done."
