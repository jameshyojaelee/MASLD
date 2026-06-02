#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=reformat_gwas
#SBATCH --output=logs/reformat_gwas_%j.out
#SBATCH --error=logs/reformat_gwas_%j.err

# Reformat our existing GWAS to finemapping pipeline input format
# Then identify lead SNPs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping

eval "$(micromamba shell hook -s bash)"
# Reformatting uses rtracklayer (only in rnaseq env)
micromamba activate rnaseq

echo "=== Step 1: Reformatting GWAS ==="
Rscript src/00_reformat_gwas.R

echo ""
echo "=== Step 2: Identifying lead SNPs ==="
bash src/01_identify_lead_snps.sh

echo ""
echo "=== All preprocessing complete ==="
