#!/bin/bash
#SBATCH --job-name=mcp_plain_nmf
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/472b_plain_nmf_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/472b_plain_nmf_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp
source "$HOME/.bashrc"; micromamba activate rapids_singlecell
python -u 472b_plain_nmf_simple.py
