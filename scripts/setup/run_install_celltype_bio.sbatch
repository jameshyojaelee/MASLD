#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --job-name=ctbio_install
#SBATCH --output=RNA-seq/logs/ctbio_install_%j.out
#SBATCH --error=RNA-seq/logs/ctbio_install_%j.err

set -eo pipefail
mkdir -p RNA-seq/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
bash scripts/setup/install_celltype_bio_env.sh
