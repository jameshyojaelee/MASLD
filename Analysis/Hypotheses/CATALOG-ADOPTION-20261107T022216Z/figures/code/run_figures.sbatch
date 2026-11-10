#!/bin/bash
#SBATCH --job-name=catalog-figures
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/figures/logs/render_%j.out
#SBATCH --error=Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/figures/logs/render_%j.err
set -euo pipefail
TASK="$PWD/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/figures"
export CATALOG_FIGURE_OUT="$TASK"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
test ! -e "$TASK/panels"
mkdir "$TASK/panels" "$TASK/source_tables" "$TASK/proofs"
"$PY" "$TASK/code/render_fig6.py" "$TASK/panels" "$TASK/source_tables"
"$RS" "$TASK/code/render_native.R"
cp -p figures/main/fig6_gene_catalog/panels/fig6a_catalog_structure.pdf figures/main/fig6_gene_catalog/panels/fig6e_portal_source_access.pdf "$TASK/panels/"
"$PY" "$TASK/code/verify_figures.py"
touch "$TASK/RENDER_COMPLETED"
