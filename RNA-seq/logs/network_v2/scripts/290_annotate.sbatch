#!/bin/bash
#SBATCH --job-name=netv2_290_annotate
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "290 | Start: $(date)"
python 290_edge_annotation_atlas.py
echo "290 | Exit: $? | End: $(date)"
