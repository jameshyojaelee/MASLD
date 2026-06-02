#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=figS09_locus_zoom
#SBATCH --output=logs/figS09_locus_zoom_%j.out
#SBATCH --error=logs/figS09_locus_zoom_%j.err

# Usage:
#   sbatch run_figS09_locus_zoom.sh                                      # default: locus_ALT_chr19_41353107 ALT
#   sbatch run_figS09_locus_zoom.sh locus_ALT_chr22_44324855 ALT         # PNPLA3 locus
#   sbatch run_figS09_locus_zoom.sh locus_GGT_chr12_111884608 GGT        # ALDH2 locus

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p logs

RSCRIPT="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
LOCUS_ID="${1:-locus_ALT_chr19_41353107}"
TRAIT="${2:-ALT}"

echo "Locus: ${LOCUS_ID} | Trait: ${TRAIT}"
echo "Start: $(date)"
echo ""

${RSCRIPT} scripts/figures/figS09_locus_zoom.R "${LOCUS_ID}" "${TRAIT}"

echo ""
echo "End: $(date)"
echo "Outputs in: GWAS/finemapping/figures/locus_zoom/"
