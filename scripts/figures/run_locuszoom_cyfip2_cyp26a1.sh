#!/bin/bash
#SBATCH --job-name=ggseqlogo
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/locuszoom_cyfip2_cyp26a1_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/locuszoom_cyfip2_cyp26a1_%j.err
#
# Regenerates the CYFIP2 + CYP26A1 4-panel Fig3-K case-study figures. Panel A
# (locus zoom) of each now writes into the consolidated Fig2 locus_zoom/ dir
# as GGT_CYFIP2.pdf / GGT_CYP26A1.pdf (2026-07-01); panels B/C/D stay in their
# existing figS_cyfip2_locus/ and figS_cyp26a1_locus/ dirs.
set +e
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE"
mkdir -p scripts/figures/logs

Rscript scripts/figures/figS_cyfip2_locus_panels.R
Rscript scripts/figures/figS_cyp26a1_locus_panels.R
echo "DONE $(date)"
