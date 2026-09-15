#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_anchors
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p2-anchors-20260909T191328Z/logs/anchors-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design AGA_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p2-anchors-20260909T191328Z PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas
[ -s "$AGA_OUT_ROOT/tables/anchor_hg38.tsv" ] || /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript 21a_anchor_liftover.R
/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python 21_anchor_panel.py
echo ANCHORS_DONE
