#!/usr/bin/env bash
#SBATCH --job-name=ag_atlas_downstream
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/logs/downstream-%j.out
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python

echo "== P1 extensions E1/E2/E3 on the full archive"
AGA_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z $PY 34_p1_extensions.py

echo "== P6a drug-target direction on the full archive"
AGA_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z $PY 60_target_direction.py

echo "== P3 allelic re-run against the FULL liver summaries (was the preview's direct-only archive)"
for c in gse281367 gse244832; do
  AGA_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z AGA_LIVER_SUMMARIES=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables/liver_summaries.tsv.gz $PY 73_atac_allelic_analysis.py $c
done
echo DOWNSTREAM_DONE
