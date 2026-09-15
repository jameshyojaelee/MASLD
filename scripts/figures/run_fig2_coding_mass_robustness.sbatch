#!/bin/bash
#SBATCH --job-name=fig2_robust
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/fig2_robust_%j.out
#SBATCH --error=logs/fig2_robust_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
P=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$P" scripts/figures/fig2_coding_mass_robustness.py \
  --project-root /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design \
  --output-dir GWAS/finemapping/results/fig2_robustness/20260907T105813Z
echo "FIG2_CODING_MASS_ROBUSTNESS_COMPLETE"
