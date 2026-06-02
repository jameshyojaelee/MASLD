#!/bin/bash
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/ncRNA_M5_epigenomic_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/ncRNA_M5_epigenomic_%j.err
# Wave1 [D] ncRNA M5: 56_ncrna_epigenomic (after M1). NOTE: uses rnaseq env (data.table/
# GenomicRanges/rtracklayer) -- NOT motifbreakr (verified: no motifbreakR/BSgenome refs).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "=== ncRNA M5 56_ncrna_epigenomic | $(date) ==="
Rscript 56_ncrna_epigenomic.R
if [ $? -ne 0 ]; then echo "FAILED: 56_ncrna_epigenomic"; exit 1; fi
echo "=== done $(date) ==="
