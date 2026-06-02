#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --job-name=ctbio_retry
#SBATCH --output=RNA-seq/logs/ctbio_retry_%j.out
#SBATCH --error=RNA-seq/logs/ctbio_retry_%j.err
set -eo pipefail
mkdir -p RNA-seq/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
bash scripts/setup/retry_celltype_bio_install.sh
