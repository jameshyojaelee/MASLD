#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=48:00:00
#SBATCH --job-name=TOAST
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/celltype_85_86_88_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/celltype_85_86_88_%j.err
# Wave1 [F]: 85 secretome chain + 86 cellstate signatures + 88 gwas-variant celltype.
# Serial (one Rscript at a time, no &). 85/88 read 80's celltype_primary_attribution.csv;
# all read -s2 dream_results_ashr.csv. NOT atlas feeders. Runs after 80.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== 85 secretome plasma chain | $(date) ==="
Rscript RNA-seq/85_secretome_plasma_chain.R
if [ $? -ne 0 ]; then echo "FAILED: 85_secretome_plasma_chain"; exit 1; fi
echo "=== 86 cellstate signatures bulk | $(date) ==="
Rscript RNA-seq/86_cellstate_signatures_bulk.R
if [ $? -ne 0 ]; then echo "FAILED: 86_cellstate_signatures_bulk"; exit 1; fi
echo "=== 88 gwas variant celltype architecture | $(date) ==="
Rscript RNA-seq/88_gwas_variant_celltype_architecture.R
if [ $? -ne 0 ]; then echo "FAILED: 88_gwas_variant_celltype_architecture"; exit 1; fi
echo "=== done $(date) ==="
