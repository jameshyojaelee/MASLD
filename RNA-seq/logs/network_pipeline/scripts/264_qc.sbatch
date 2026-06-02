#!/bin/bash
#SBATCH --job-name=net_264_qc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "264 network_qc | Start: $(date)"
python 264_network_qc.py
echo "264 network_qc | Exit: $? | End: $(date)"
