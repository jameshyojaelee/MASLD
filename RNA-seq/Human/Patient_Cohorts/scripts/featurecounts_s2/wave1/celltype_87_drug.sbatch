#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=TOAST
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/celltype_87_drug_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/celltype_87_drug_%j.err
# Wave1 [F]: 87 celltype drug-target enrichment. Reads 80's attribution + -s2 dream_ashr +
# drug-v3 outputs (opentargets_known_drugs.csv, top50). Runs after 80 AND drug-v3. NOT a feeder.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== 87 celltype drug target enrichment | $(date) ==="
Rscript RNA-seq/87_celltype_drug_target_enrichment.R
if [ $? -ne 0 ]; then echo "FAILED: 87_celltype_drug_target_enrichment"; exit 1; fi
echo "=== done $(date) ==="
