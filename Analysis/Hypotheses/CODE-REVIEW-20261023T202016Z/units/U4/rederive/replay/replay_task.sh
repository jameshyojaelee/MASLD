#!/bin/bash
# Replay one sampled (study, chr, gene) pair with the unmodified producer
# 06_susie_coloc.R, in the environment of 09b_coloc_eligibility_batch.sbatch.
# Usage: replay_task.sh <task index in pairs.tsv>
set -eo pipefail
K=$1
RP=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/replay
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
LDREF=${FM}/data/ld_ref
ROW=$(awk -F'\t' -v k="$K" 'NR>1 && $1==k' ${RP}/pairs.tsv)
GWAS=$(echo "$ROW" | cut -f3); CHR=$(echo "$ROW" | cut -f4); ENSG=$(echo "$ROW" | cut -f5)
TD=${RP}/out/task_${K}
mkdir -p ${TD}
printf "gwas_name\tchr\tensembl\n%s\t%s\t%s\n" "$GWAS" "$CHR" "$ENSG" > ${TD}/filter.tsv
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AFR_LD_DIR="${LDREF}/1kg_afr_gram" AMR_LD_DIR="${LDREF}/1kg_amr_gram"
export EAS_LD_DIR="${LDREF}/1kg_eas_gram" SAS_LD_DIR="${LDREF}/1kg_sas_gram"
export EQTL_SUSIE_DIR="${FM}/results/eqtl_susie_polyfun"
export COLOC_OUT_ROOT="${TD}/coloc"
export COLOC_GENE_FILTER="${TD}/filter.tsv"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "[$(date -u +%FT%TZ)] task ${K}: ${GWAS} chr${CHR} ${ENSG}"
sha256sum ${FM}/src/06_susie_coloc.R ${FM}/src/coloc_signal_eligibility.R ${FM}/src/perf/get_ld_per_locus_cached.R ${FM}/src/finemapping_functions.R > ${TD}/code_sha256.tsv
START=$(date +%s)
cd ${FM}/src
GWAS_NAME="${GWAS}" CHR_FILTER="${CHR}" /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript ${FM}/src/06_susie_coloc.R > ${TD}/run.log 2>&1
echo "[$(date -u +%FT%TZ)] task ${K} done in $(( $(date +%s) - START )) s"
