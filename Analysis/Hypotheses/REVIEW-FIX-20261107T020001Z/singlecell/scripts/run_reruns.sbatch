#!/bin/bash
#SBATCH --job-name=sc-lineage-fix
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/singlecell/logs/rerun_%j.out
#SBATCH --error=Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/singlecell/logs/rerun_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
TASK="$PWD/Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/singlecell"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export R_DATATABLE_NUM_THREADS=1
test ! -e "$TASK/results"
mkdir "$TASK/results"
"$RS" "$TASK/scripts/check_lineage_cells.R"
"$RS" -e 'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("MASLD_PROJECT_ROOT"),"Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/singlecell/results/sessionInfo.txt"))'
export SINGLECELL_BINARY_DE_OUT="$TASK/results/binary"
"$RS" "$TASK/scripts/pseudobulk_de.R"
export SINGLECELL_BINARY_DE_OUT="$TASK/results/strict_control"
export SINGLECELL_STRICT_CONTROL=1
"$RS" "$TASK/scripts/pseudobulk_de.R"
unset SINGLECELL_STRICT_CONTROL
export SINGLECELL_STAGE_DE_OUT="$TASK/results/coarse"
"$RS" "$TASK/scripts/run_pseudobulk_de_coarse_stages.R"
export SINGLECELL_BINARY_DE_OUT="$TASK/results/binary"
export SINGLECELL_ALLCELL_OUT="$TASK/results/allcell"
"$RS" "$TASK/scripts/314_allcell_pseudobulk_de.R"
touch "$TASK/RERUNS_COMPLETED"
