#!/bin/bash
#SBATCH --job-name=net_294_portal_rebuild
#SBATCH --partition=bigmem
#SBATCH --mem=240G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=logs/net_294_portal_rebuild_%j.out
#SBATCH --error=logs/net_294_portal_rebuild_%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

echo "294 portal rebuild | Start: $(date)"
python 294_export_portal_v2.py
echo "294 portal rebuild | Exit: $? | End: $(date)"
