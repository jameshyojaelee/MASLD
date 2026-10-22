#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN5/refit/%j.out
# AN5 task 2: GWAS-side SuSiE diagnostics + refits for the 14 PP.H4 = 1.0 genes (18 rows) and 3 comparator rows.
# Runs a COPY of GWAS/finemapping/src/06_susie_coloc.R with a diagnostic block inserted (refit/an5_diag_block.R),
# same environment as src/perf/09b_coloc_eligibility_batch.sbatch (polyfun eQTL fits, Gram non-EUR LD). Seeds 42 as in the source.
set -eo pipefail
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
FM=${MASLD_PROJECT_ROOT}/GWAS/finemapping
LDREF=${FM}/data/ld_ref
export AFR_LD_DIR=${LDREF}/1kg_afr_gram AMR_LD_DIR=${LDREF}/1kg_amr_gram EAS_LD_DIR=${LDREF}/1kg_eas_gram SAS_LD_DIR=${LDREF}/1kg_sas_gram
export EQTL_SUSIE_DIR=${FM}/results/eqtl_susie_polyfun
export COLOC_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN5/refit/run_20261021T125037Z/coloc
export COLOC_GENE_FILTER=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN5/refit/run_20261021T125037Z/plan_rows.tsv
export AN5_DIAG_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN5/refit/run_20261021T125037Z/diag
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript -e 'library(coloc); library(susieR); library(data.table); sessionInfo()' > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN5/refit/run_20261021T125037Z/env/sessionInfo.txt 2>&1
sha256sum ${FM}/src/06_susie_coloc.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN5/refit/06_susie_coloc_AN5diag.R ${FM}/src/coloc_signal_eligibility.R > /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN5/refit/run_20261021T125037Z/env/script_sha256.txt
tail -n +2 /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN5/refit/run_20261021T125037Z/plan_rows.tsv | cut -f1,2 | sort -u | while IFS=$'\t' read -r GWAS CHR; do
  echo "[$(date -u +%FT%TZ)] ${GWAS} chr${CHR}"
  GWAS_NAME="${GWAS}" CHR_FILTER="${CHR}" Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN5/refit/06_susie_coloc_AN5diag.R || echo "TASK FAILED ${GWAS} chr${CHR}"
done
echo "[$(date -u +%FT%TZ)] done"
