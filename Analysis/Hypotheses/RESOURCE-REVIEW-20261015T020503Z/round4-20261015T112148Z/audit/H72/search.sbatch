#!/bin/bash
#SBATCH --job-name=9723
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H72/9723_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H72/9723_%j.err
set -uo pipefail
R=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
A=$R/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H72
export TMPDIR=$A/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
cd $R
PAT='PanUKBB_F_|PanUKBB_M_|female\.varorder|male\.varorder|30620_raw|30650_raw|30730_raw|z_diff|zdiff|beta_F|sex_het|sex-heterogen|sexhet|sex_interaction|GxSex'
# targeted outcome search; environments, package trees and the H72 audit dir are pruned
for D in Analysis/Hypotheses RNA-seq GWAS/MR_Data GWAS/finemapping/scripts GWAS/finemapping/results scripts data/public_r4; do
  find $D \( -name '.mamba' -o -name 'site-packages' -o -name 'conda-meta' -o -name 'environments' -o -name 'env' -o -path '*/audit/H72' \) -prune -o \
    -type f -size -20M \( -name '*.py' -o -name '*.R' -o -name '*.sh' -o -name '*.tsv' -o -name '*.csv' -o -name '*.txt' -o -name '*.md' -o -name '*.json' -o -name '*.out' -o -name '*.log' \) -print0 2>/dev/null \
  | xargs -0 -P 4 grep -IlE "$PAT" 2>/dev/null
done | sort > $A/outcome_search_targeted.txt
echo "search done $(wc -l < $A/outcome_search_targeted.txt) files" > $A/outcome_search_targeted.done
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze > $A/pip_freeze.txt 2>/dev/null
sha256sum $A/audit_h72.py $A/audit.sbatch $A/audit_power_281.py $A/audit_power.sbatch $A/search.sbatch > $A/audit_script_sha256.txt
