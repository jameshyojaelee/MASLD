#!/usr/bin/env bash
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=${ROOT}/GWAS/finemapping/src/seqfunc_v2/chrombpnet
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2
SMOKE_JOB=${1:?usage: submit_training_v2.sh SMOKE_JOB_ID}

# Fold0 is the technical/internal gate. Remaining folds cannot start unless its
# full package-native held-out metrics meet every frozen threshold.
BIAS0=$(sbatch --parsable --array=0 --dependency=afterok:${SMOKE_JOB} "${SRC}/07_train_bias.sbatch")
MODEL0=$(sbatch --parsable --array=0 --dependency=afterok:${BIAS0} "${SRC}/08_train_model.sbatch")
QC0=$(sbatch --parsable --array=0 --dependency=afterok:${MODEL0} "${SRC}/09_internal_qc.sbatch")
BIASREST=$(sbatch --parsable --array=1-4%4 --dependency=afterok:${QC0} "${SRC}/07_train_bias.sbatch")
MODELREST=$(sbatch --parsable --array=1-4%4 --dependency=afterok:${BIASREST} "${SRC}/08_train_model.sbatch")
QCREST=$(sbatch --parsable --array=1-4%4 --dependency=afterok:${MODELREST} "${SRC}/09_internal_qc.sbatch")
printf 'stage\tjob_id\tdependency\nfold0_bias\t%s\tafterok:%s\nfold0_model\t%s\tafterok:%s\nfold0_qc\t%s\tafterok:%s\nfold1_4_bias\t%s\tafterok:%s\nfold1_4_model\t%s\tafterok:%s\nfold1_4_qc\t%s\tafterok:%s\n' \
  "${BIAS0}" "${SMOKE_JOB}" "${MODEL0}" "${BIAS0}" "${QC0}" "${MODEL0}" \
  "${BIASREST}" "${QC0}" "${MODELREST}" "${BIASREST}" "${QCREST}" "${MODELREST}" \
  | tee "${OUT}/submitted_training_jobs.tsv"
