#!/bin/bash
#SBATCH --job-name=NMF
#SBATCH --partition=cpu,io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/208_rerun_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/208_rerun_%j.err
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "=== 208 subtype_coloc rerun (sex file provisioned) $(date) ==="
Rscript 208_subtype_coloc.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 208"; exit 1; fi
echo "Finished: $(date)"
