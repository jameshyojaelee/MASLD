#!/bin/bash
#SBATCH --job-name=net_265_assemble
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00


eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "265 265_assemble_multiplex_graph.py | Start: $(date)"
python 265_assemble_multiplex_graph.py
echo "265 265_assemble_multiplex_graph.py | Exit: $? | End: $(date)"
