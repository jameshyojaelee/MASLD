#!/usr/bin/env bash
#SBATCH --job-name=agp_f2-haplotype-report
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f2-haplotype-report-%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export F2_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program/f2-haplotype-20260914T230433Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program
"$PY" f2_08_write_results.py
"$PY" f2_04_manifest.py
"$PY" -m pip freeze > "$F2_OUT_ROOT/pip_freeze.txt" 2>/dev/null || true
cp f2_recipe.py f2_01_build_null_panel.py f2_01_liftover.R f2_02_score.py f2_03_analyse.py f2_04_manifest.py \
   f2_06_arm_frequencies.py f2_07_deposit_readout_check.py f2_08_write_results.py \
   f2_01_panel.sbatch f2_02_score.sbatch f2_02_score_ext2500.sbatch f2_05_gnmt_ld.sbatch \
   f2_06_arm_frequencies.sbatch f2_07_readout_check.sbatch f2_09_final.sbatch "$F2_OUT_ROOT/code_snapshot/"
echo F2_REPORT_DONE
