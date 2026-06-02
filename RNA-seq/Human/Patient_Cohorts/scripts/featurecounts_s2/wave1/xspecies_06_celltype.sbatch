#!/bin/bash
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/xspecies_06_celltype_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/xspecies_06_celltype_%j.err
# Wave1 [E] TERMINAL: 06 cell-type concordance (human scRNA pseudobulk vs mouse bulk dream_pooled).
# Independent of 01-05; NOT read by the 05 unified atlas. Reads mouse meta_analysis dream_pooled.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
S=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== xspecies 06 cell-type concordance | $(date) ==="
Rscript "$S/06_celltype_concordance.R"
if [ $? -ne 0 ]; then echo "FAILED: 06_celltype_concordance"; exit 1; fi
echo "=== done $(date) ==="
