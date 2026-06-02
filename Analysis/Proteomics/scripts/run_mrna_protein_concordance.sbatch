#!/bin/bash
#SBATCH --job-name=mrna_prot_conc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/mrna_protein_concordance_%j.log
#SBATCH --error=logs/mrna_protein_concordance_%j.log

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Proteomics
mkdir -p logs

echo "=== mRNA-protein concordance (Gap #10) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Start: $(date)"

eval "$(micromamba shell hook -s bash)"
set +u  # conda activation references unbound vars (ADDR2LINE) — would abort under set -u
micromamba activate rnaseq
set -u

Rscript scripts/mrna_protein_concordance.R

echo "End: $(date)"
