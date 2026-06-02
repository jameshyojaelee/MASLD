#!/bin/bash
#SBATCH --job-name=install_deps
#SBATCH --output=logs/install_deps_%j.out
#SBATCH --error=logs/install_deps_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

source /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/scripts/environment_setup.sh


echo "Installing missing dependencies..."
micromamba install -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/.mamba/pathway_analysis -y -c conda-forge -c bioconda \
    bioconductor-gsva \
    r-data.table


echo "Installation complete."
