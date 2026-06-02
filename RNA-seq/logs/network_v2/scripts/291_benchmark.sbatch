#!/bin/bash
#SBATCH --job-name=netv2_291_benchmark
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "291 | Start: $(date)"
Rscript 291_edge_quality_report.R
echo "291 | Exit: $? | End: $(date)"
