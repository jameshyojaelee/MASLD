#!/bin/bash
#SBATCH --job-name=net_276_lls
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
python 276_lls_calibration.py
