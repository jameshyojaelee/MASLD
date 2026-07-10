#!/bin/bash
#SBATCH --job-name=matplotlib
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=72G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/hep_hotspot_composite_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/hep_hotspot_composite_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/figures/logs
micromamba run -n spatial python scripts/figures/hep_hotspot_composite.py
