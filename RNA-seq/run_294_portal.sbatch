#!/bin/bash
#SBATCH --job-name=294_portal_export
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/294_portal_export_%j.out
#SBATCH --error=logs/294_portal_export_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
micromamba run -n spatial python 294_export_portal_v2.py
