#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=4 --mem=16G --time=6:00:00
#SBATCH --job-name=bioc_only --output=RNA-seq/logs/bioc_only_%j.out --error=RNA-seq/logs/bioc_only_%j.err
set -eo pipefail; mkdir -p RNA-seq/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
bash scripts/setup/bioc_only_install.sh
