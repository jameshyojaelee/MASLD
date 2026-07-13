#!/usr/bin/env bash
#SBATCH --job-name=MACS3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-2%3
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/macs3_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/macs3_%A_%a.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2
BLACKLIST=/gpfs/commons/home/jameslee/reference_genome/blacklists/hg38-blacklist.v2.bed
SAMPLES=(pooled pseudorep1 pseudorep2)
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
OUT=${BASE}/peaks/${SAMPLE}
TAG=${BASE}/fragments/hepatocyte.${SAMPLE}.tn5.tagAlign.gz
mkdir -p "${OUT}"
test -s "${TAG}"

micromamba run -p "${ENV}" macs3 callpeak \
  -t "${TAG}" \
  -f BED \
  -g hs \
  -n hepatocyte.${SAMPLE} \
  --outdir "${OUT}" \
  --nomodel \
  --shift -100 \
  --extsize 200 \
  --keep-dup all \
  -q 0.05

RAW=${OUT}/hepatocyte.${SAMPLE}_peaks.narrowPeak
FILTERED=${OUT}/hepatocyte.${SAMPLE}.autosomal.blacklist_filtered.narrowPeak
test -s "${RAW}"
# Preserve all ten genuine MACS3 narrowPeak columns, including signalValue,
# qValue, and summit offset. Exclude any blacklist overlap.
micromamba run -p "${ENV}" bedtools intersect -v -a "${RAW}" -b "${BLACKLIST}" \
  | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ && $2 >= 0 && $3 > $2 && NF >= 10' \
  > "${FILTERED}.tmp"
test -s "${FILTERED}.tmp"
mv "${FILTERED}.tmp" "${FILTERED}"
awk 'BEGIN{OFS="\t"} {print $0, $2+$10, $2+$10+1}' "${FILTERED}" \
  > "${OUT}/hepatocyte.${SAMPLE}.peak_summits.audit.tsv"
printf 'sample\traw_peaks\tfiltered_peaks\n%s\t%s\t%s\n' \
  "${SAMPLE}" "$(wc -l < "${RAW}")" "$(wc -l < "${FILTERED}")" \
  > "${OUT}/peak_call_qc.tsv"

