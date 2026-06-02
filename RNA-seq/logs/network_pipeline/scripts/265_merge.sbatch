#!/bin/bash
#SBATCH --job-name=net_265_merge
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "265 graph_merge | Start: $(date)"
python 265_graph_merge.py
echo "265 graph_merge | Exit: $? | End: $(date)"
