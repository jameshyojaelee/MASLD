#!/bin/bash
#SBATCH --job-name=net_290_annot_rebuild
#SBATCH --partition=bigmem
#SBATCH --mem=256G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=logs/net_290_annot_rebuild_%j.out
#SBATCH --error=logs/net_290_annot_rebuild_%j.err

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

echo "290 atlas rebuild | Start: $(date)"
python 290_edge_annotation_atlas.py
echo "290 atlas rebuild | Exit: $? | End: $(date)"
