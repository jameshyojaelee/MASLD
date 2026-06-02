#!/bin/bash
#SBATCH --job-name=net_270_report
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "270 pipeline_report | Start: $(date)"
python 270_pipeline_report.py
echo "270 pipeline_report | Exit: $? | End: $(date)"
