#!/bin/bash
#SBATCH --job-name=B4_step4_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/logs/B4_step4_coloc_%j.out
#SBATCH --error=Analysis/SingleCell/logs/B4_step4_coloc_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/logs RNA-seq/results/stratified_causal

# rnaseq env's binutils activate hook references unbound ADDR2LINE; relax `-u` around activate.
set +u
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

Rscript Analysis/SingleCell/scripts/353_celltype_sex_coloc.R

echo "Step 4 done: $(date)"
