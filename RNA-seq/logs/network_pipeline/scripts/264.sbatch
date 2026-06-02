#!/bin/bash
#SBATCH --job-name=net_264_diagnost
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00


eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "264 264_bayesian_diagnostics.py | Start: $(date)"
python 264_bayesian_diagnostics.py
echo "264 264_bayesian_diagnostics.py | Exit: $? | End: $(date)"
