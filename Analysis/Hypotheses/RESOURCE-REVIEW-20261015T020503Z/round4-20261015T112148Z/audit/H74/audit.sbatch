#!/bin/bash
#SBATCH --job-name=9741
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H74/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H74/logs/%j.err
set -euo pipefail
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H74
RR=Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z
OUT=${A}/run-$(date -u +%Y%m%dT%H%M%SZ)
test ! -e "${OUT}"; mkdir -p "${OUT}"
export TMPDIR=${A}/tmp XDG_CACHE_HOME=${A}/tmp PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
for f in GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv \
         data/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf \
         GWAS/finemapping/config/gwas_trait_tier.tsv \
         ${RR}/gate/NG13/count.py ${RR}/gate/NG13/count.sbatch ${RR}/gate/NG13/count_22251943.log \
         ${RR}/gate/NG13/locus_structure.tsv ${RR}/write/H74/count_groups.sh ${RR}/write/H74/count_groups.out; do
  printf "%s\t%s\t%s\n" "$f" "$(stat -c %s $f)" "$(sha256sum $f | cut -d' ' -f1)"
done > ${OUT}/input_sha256_recomputed.tsv
bash ${RR}/write/H74/count_groups.sh > ${OUT}/count_groups_rerun.out
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python ${A}/audit_counts.py ${OUT} > ${OUT}/audit_counts.log 2>&1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze > ${OUT}/pip_freeze.txt
sha256sum ${A}/audit_counts.py ${A}/audit.sbatch > ${OUT}/script_sha256.tsv
