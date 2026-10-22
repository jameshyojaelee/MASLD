#!/bin/bash
#SBATCH --job-name=9721
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H72/9721_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H72/9721_%j.err
set -uo pipefail
R=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
A=$R/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H72
export TMPDIR=$A/tmp XDG_CACHE_HOME=$A/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
cd $R
# outcome search: any on-disk file (<20 MB, text) that compares female and male enzyme effects
( find Analysis RNA-seq GWAS/MR_Data GWAS/finemapping/scripts GWAS/finemapping/results scripts -type f -size -20M \
    \( -name '*.py' -o -name '*.R' -o -name '*.sh' -o -name '*.tsv' -o -name '*.csv' -o -name '*.txt' -o -name '*.md' -o -name '*.json' -o -name '*.out' -o -name '*.log' \) \
    -not -path '*/audit/H72/*' -print0 2>/dev/null \
  | xargs -0 grep -IlE 'PanUKBB_F_|PanUKBB_M_|female\.varorder|male\.varorder|30620_raw|30650_raw|30730_raw|z_diff|zdiff|beta_F|sex_het|sex-heterogen|sexhet|sex_interaction|GxSex|sex_strat' 2>/dev/null \
  > $A/outcome_search_files.txt ) &
$PY $A/audit_h72.py
wait
$PY -m pip freeze > $A/pip_freeze.txt 2>/dev/null
sha256sum $A/audit_h72.py $A/audit.sbatch > $A/audit_script_sha256.txt
