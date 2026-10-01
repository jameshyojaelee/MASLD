#!/bin/bash
#SBATCH --job-name=masld_query_figure
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_query_figure_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_query_figure_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python \
 "$PROJ/scripts/analysis/alphagenome_campaign/coordinator_query_figure.py" \
 --db "$BASE/catalog-currin-21773574/regulatory.sqlite" \
 --example "$BASE/catalog-currin-21773574/example_currin_variant.json" \
 --check "$BASE/review/catalog-currin-21773580/checks.json" \
 --out "$PROJ/figures/candidates/regulatory-query-week1-20260915-v1"
