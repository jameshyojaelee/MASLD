#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --job-name=B1ag
#SBATCH --output=RNA-seq/logs/B1ag_%j.out
#SBATCH --error=RNA-seq/logs/B1ag_%j.err

# B1 F/M aggregator: pools PanUKBB_F_* + PanUKBB_M_* SuSiE-COLOC chrs,
# produces RNA-seq/results/stratified_causal/sex_stratified_coloc_panukbb.csv
# with sex-specificity scores. Tolerates missing strata (M_AST, M_GGT cancelled).

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/logs RNA-seq/results/stratified_causal

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "Start: $(date)"
Rscript RNA-seq/216a_aggregate_sex_strat_coloc.R
echo "End: $(date)"
