#!/bin/bash
#SBATCH --job-name=u4_abf_sdy
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/logs/%j.out
# ABF arm under estimated sdY (eQTL 0.714; GWAS chr22 regression estimates) for
# 6,716 tier-1/2 pairs with released PP.H4.abf in [0.1, 0.95]; 751 (study, chr) tasks.
# Task 1 (chr1, largest eQTL table) is timed first; parallelism is set from its peak RSS.
set -eo pipefail
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/abf_sdy
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sha256sum ${A}/06_susie_coloc_sdy.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/src/06_susie_coloc.R ${A}/tasks.tsv > ${A}/code_sha256_${SLURM_JOB_ID}.tsv
/usr/bin/time -v -o ${A}/time_first_task.txt bash ${A}/run_line.sh 2
PEAK_KB=$(awk -F': ' '/Maximum resident/{print $2}' ${A}/time_first_task.txt)
P=$(( 58000000 / (PEAK_KB * 5 / 4) )); (( P > 15 )) && P=15; (( P < 1 )) && P=1
echo "first task peak ${PEAK_KB} kB -> ${P} parallel processes"
seq 3 $(wc -l < ${A}/tasks.tsv) | xargs -P ${P} -I{} bash ${A}/run_line.sh {}
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla ${A}/compare.R > ${A}/compare.log 2>&1
echo "[$(date -u +%FT%TZ)] all done"
