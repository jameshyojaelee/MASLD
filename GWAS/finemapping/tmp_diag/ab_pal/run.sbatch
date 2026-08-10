#!/bin/bash
#SBATCH --job-name=susie_coloc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/ab_pal/%x_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/ab_pal/%x_%A_%a.err
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate finemapping
set -u
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
cd "$FM"
# array index -> (study, chr, arm). ARM A = palindromes RETAINED (new code),
# ARM B = COLOC_DROP_PALINDROMIC=1 (old behaviour). Same eQTL fits, same LD,
# so the ONLY difference is the palindrome change.
STUDIES=(MVP_NAFLD_EUR MVP_ALT_EUR MVP_ALT_AFR)
i=$SLURM_ARRAY_TASK_ID
s=$(( (i / 2) % 3 )); arm=$(( i % 2 ))
STUDY=${STUDIES[$s]}
export GWAS_NAME=$STUDY CHR_FILTER=22
if [ "$arm" -eq 1 ]; then
  export COLOC_DROP_PALINDROMIC=1
  export COLOC_OUT_SUFFIX=_abtest_drop
else
  export COLOC_OUT_SUFFIX=_abtest_keep
fi
echo "=== $STUDY chr22 arm=$arm suffix=$COLOC_OUT_SUFFIX ==="
Rscript src/06_susie_coloc.R
