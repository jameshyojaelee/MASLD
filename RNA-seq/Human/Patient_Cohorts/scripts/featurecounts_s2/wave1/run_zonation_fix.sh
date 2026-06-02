#!/bin/bash
#SBATCH --job-name=fgsea
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/zonation_fix_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/zonation_fix_%j.err
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "=== 43_zonation (STAR -s2, schema-robust nafl_nash fix) $(date) ==="
Rscript 43_zonation_classification.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 43_zonation"; exit 1; fi
echo "Finished: $(date)"
