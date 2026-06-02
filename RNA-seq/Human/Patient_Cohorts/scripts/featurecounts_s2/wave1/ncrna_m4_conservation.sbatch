#!/bin/bash
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/ncRNA_M4_conservation_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/ncRNA_M4_conservation_%j.err
# Wave1 [D] ncRNA M4: 55_ncrna_conservation (after M1). Synteny conservation.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "=== ncRNA M4 55_ncrna_conservation | $(date) ==="
Rscript 55_ncrna_conservation.R
if [ $? -ne 0 ]; then echo "FAILED: 55_ncrna_conservation"; exit 1; fi
echo "=== done $(date) ==="
