#!/bin/bash
#SBATCH --job-name=B1_kallisto_agg
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=logs/agg_%j.out
#SBATCH --error=logs/agg_%j.err

set -eo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

WT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
cd "$WT"
Rscript RNA-seq/scripts/kallisto/03_aggregate_tximport.R
