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
PREPARED=${CANDIDATE}/genetics/prepared
PLAN=${PREPARED}/replay_batches/batch_${BATCH}.tsv
RUNNER=${PREPARED}/replay_runner/06_susie_coloc_replay.R
HELPER=${PREPARED}/replay_runner/genetics_export_helpers.R
GATE=${PREPARED}/PROMOTED_UPSTREAM_GATE.tsv
PENDING=${CANDIDATE}/genetics/.batch_${BATCH}.pending.${SLURM_JOB_ID}
FINAL=${CANDIDATE}/genetics/replay_execution/batch_${BATCH}
RUN_ROOT=${PENDING}/runs
EXPORT_ROOT=${PENDING}/exports
SCRATCH=$(mktemp -d /scratch/${USER}/atac_v3_replay_${BATCH}_${SLURM_JOB_ID}_XXXXXX)
trap 'rm -rf -- "${SCRATCH}"' EXIT
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
eval "$(${MAMBA} shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u
python "${PACKAGE}/21_verify_replay_batch_inputs.py" --batch "${BATCH}"
if [[ -e "${PENDING}" || -e "${FINAL}" ]]; then
  echo "Refusing to overwrite replay batch state" >&2
  exit 1
fi
mkdir -p "${PENDING}" "$(dirname "${FINAL}")"
GWAS_REGISTRY=$(awk -F '\t' '
  NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
  NR==2 {print $(h["gwas_registry_path"])}
' "${GATE}")
EQTL_DIR=$(awk -F '\t' '
  NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
  NR==2 {print $(h["eqtl_susie_dir"])}
' "${GATE}")
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
  ATAC_V3_REPLAY_HELPER="${HELPER}" \
  ATAC_V3_GWAS_REGISTRY="${GWAS_REGISTRY}" \
  EQTL_SUSIE_DIR="${EQTL_DIR}" \
  Rscript "${RUNNER}" "${STUDY}" "${CHR}"
done
python "${PACKAGE}/22_validate_replay_batch.py" \
  --batch "${BATCH}" --pending "${PENDING}"
mv -- "${PENDING}" "${FINAL}"
