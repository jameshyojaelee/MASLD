#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --array=0-4%5
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/foldprep_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/foldprep_%A_%a.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
export PATH=${ENV}/bin:${PATH}
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2
FASTA=/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa
PEAKS=${BASE}/hepatocyte.idr_pooled_summit.narrowPeak
CBP=${ENV}/bin/chrombpnet
FOLD=${SLURM_ARRAY_TASK_ID}
FOLD_DIR=${BASE}/fold${FOLD}
mkdir -p "${FOLD_DIR}"
test -s "${PEAKS}"

TEST=$(awk -F '\t' -v f="${FOLD}" 'NR>1 && $1==f {print $4}' "${BASE}/fold_manifest.tsv" | tr ',' ' ')
VALID=$(awk -F '\t' -v f="${FOLD}" 'NR>1 && $1==f {print $7}' "${BASE}/fold_manifest.tsv" | tr ',' ' ')
test -n "${TEST}"
test -n "${VALID}"

# shellcheck disable=SC2086
"${CBP}" prep splits \
  -op "${FOLD_DIR}/fold${FOLD}" \
  -c "${BASE}/GRCh38.autosomes.chrom.sizes" \
  -tcr ${TEST} \
  -vcr ${VALID}
test -s "${FOLD_DIR}/fold${FOLD}.json"

"${CBP}" prep nonpeaks \
  -g "${FASTA}" \
  -o "${FOLD_DIR}/hepatocyte.nonpeaks" \
  -p "${PEAKS}" \
  -c "${BASE}/GRCh38.autosomes.chrom.sizes" \
  -fl "${FOLD_DIR}/fold${FOLD}.json" \
  -s 42
test -s "${FOLD_DIR}/hepatocyte.nonpeaks_negatives.bed"

printf 'fold\ttest_chromosomes\tvalidation_chromosomes\tpeaks\tnegatives\n%s\t%s\t%s\t%s\t%s\n' \
  "${FOLD}" "${TEST// /,}" "${VALID// /,}" "$(wc -l < "${PEAKS}")" \
  "$(wc -l < "${FOLD_DIR}/hepatocyte.nonpeaks_negatives.bed")" \
  > "${FOLD_DIR}/prep_qc.tsv"
