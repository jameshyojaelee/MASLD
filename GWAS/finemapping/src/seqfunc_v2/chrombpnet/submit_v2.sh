#!/usr/bin/env bash
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=${ROOT}/GWAS/finemapping/src/seqfunc_v2/chrombpnet
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2
mkdir -p "${OUT}/logs"
export MASLD_PROJECT_ROOT=${ROOT}

python3 "${SRC}/00_preflight.py"
EXTRACT=$(sbatch --parsable "${SRC}/01_extract_fragments.sbatch")
MACS3=$(sbatch --parsable --dependency=afterok:${EXTRACT} "${SRC}/02_macs3.sbatch")
IDR=$(sbatch --parsable --dependency=afterok:${MACS3} "${SRC}/03_idr_filter.sbatch")
FOLDS=$(sbatch --parsable --dependency=afterok:${IDR} "${SRC}/04_prepare_folds.sbatch")
SMOKE=$(sbatch --parsable --dependency=afterok:${FOLDS} "${SRC}/05_fold0_bias_smoke.sbatch")
printf 'stage\tjob_id\tdependency\nextract\t%s\t\nmacs3\t%s\tafterok:%s\nidr\t%s\tafterok:%s\nfoldprep\t%s\tafterok:%s\nbias_smoke\t%s\tafterok:%s\n' \
  "${EXTRACT}" "${MACS3}" "${EXTRACT}" "${IDR}" "${MACS3}" "${FOLDS}" "${IDR}" "${SMOKE}" "${FOLDS}" \
  | tee "${OUT}/submitted_jobs.tsv"

