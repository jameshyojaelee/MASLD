#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=mesusietopldfull
#SBATCH --output=GWAS/finemapping/logs/mesusietopldfull_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/mesusietopldfull_%A_%a.err
#
# Run MESuSiE with BOTH EUR and EAS arms on TOP-LD (parallel to SuSiEX 7c
# "topld_full"). EAS arm uses topld_eas — MESuSiE only needs .ld + .bim per
# block, which topld_eas already has (no chr-level PLINK build required, unlike
# SuSiEX which needs chr-level .bed/.fam for EAS).
#
# Usage (SLURM array over all shared loci):
#   sbatch --array=1-140 src/run_mesusie_topld_full.sh

set -eo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

LOCUS_ROW=${SLURM_ARRAY_TASK_ID:-1}
export LD_PANEL=topld
export MESUSIE_RESULTS_SUFFIX="_topld_full"
# Leave EAS_LD_DIR UNSET so dispatch resolves to topld_eas via finemapping_functions.R
unset EAS_LD_DIR AFR_LD_DIR SAS_LD_DIR UKBB_LD_DIR

echo "[mesusietopldfull] LOCUS_ROW=${LOCUS_ROW}  LD_PANEL=${LD_PANEL}  EAS_LD_DIR=<unset> (→topld_eas)  start $(date)"
Rscript src/10b_run_mesusie.R
echo "[mesusietopldfull] LOCUS_ROW=${LOCUS_ROW}  done  $(date)"
