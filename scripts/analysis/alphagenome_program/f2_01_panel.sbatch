#!/usr/bin/env bash
#SBATCH --job-name=agp_f2-haplotype-panel
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f2-haplotype-panel-%j.out
set -euo pipefail
module load PLINK/1.9b_6.21-x86_64
export PYTHONNOUSERSITE=1
export F2_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program/f2-haplotype-20260914T230433Z
export F2_SCRATCH="$F2_OUT_ROOT/raw/null_panel_build"
mkdir -p "$F2_SCRATCH"
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python f2_01_build_null_panel.py 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 19 20 22
rm -rf "$F2_SCRATCH"
echo F2_PANEL_DONE
