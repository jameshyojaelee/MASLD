#!/bin/bash
#SBATCH --job-name=atacreplay
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_replay_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_replay_%j.err

set -euo pipefail
if [[ $# -ne 1 || ! $1 =~ ^[1-5]$ ]]; then
  echo "Usage: sbatch 08_genetics_batch.sbatch 1|2|3|4|5" >&2
  exit 2
fi
BATCH=$1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
CANDIDATE=${ROOT}/Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1
PLAN=${CANDIDATE}/genetics/replay_batches/batch_${BATCH}.tsv
RUNNER=${CANDIDATE}/genetics/replay_runner/06_susie_coloc_replay.R
RUN_ROOT=${CANDIDATE}/genetics/replay_runs/batch_${BATCH}
EXPORT_ROOT=${CANDIDATE}/genetics/replay_exports
SCRATCH=$(mktemp -d /scratch/${USER}/atac_v3_replay_${BATCH}_${SLURM_JOB_ID}_XXXXXX)
trap 'rm -rf -- "${SCRATCH}"' EXIT
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
eval "$(${MAMBA} shell hook --shell bash)"
micromamba activate rnaseq
tail -n +2 "${PLAN}" | cut -f3,6 | sort -u | while IFS=$'\t' read -r STUDY CHR; do
  FILTER=${SCRATCH}/${STUDY}.chr${CHR}.tsv
  awk -F '\t' -v study="${STUDY}" -v chr="${CHR}" '
    BEGIN {OFS="\t"}
    NR==1 {for(i=1;i<=NF;i++) h[$i]=i; print "ensembl"; next}
    $(h["gwas_name"])==study && $(h["chr"])==chr {print $(h["ensembl"])}
  ' "${PLAN}" > "${FILTER}"
  ATAC_V3_REPLAY_GENE_FILTER="${FILTER}" \
  ATAC_V3_REPLAY_RUN_ROOT="${RUN_ROOT}" \
  ATAC_V3_REPLAY_EXPORT_ROOT="${EXPORT_ROOT}" \
  Rscript "${RUNNER}" "${STUDY}" "${CHR}"
done
