#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name=na_fix_test

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping

echo "=== Test 1: FinnGen_NAFLD chr4 locus (previously 0 variants after NA removal) ==="
Rscript src/03_run_fm_per_locus.R FinnGen_NAFLD EUR 4.88231865 438857 4614 0.5 EUR

echo ""
echo "=== Test 2: UKBB_AST chr1 locus (variant ordering bug) ==="
Rscript src/03_run_fm_per_locus.R UKBB_AST EUR 1.150295131 343850 0 0.5 EUR
