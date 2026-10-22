#!/bin/bash
#SBATCH --job-name=comp_774
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/bulk_refit_followup/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/bulk_refit_followup/%x_%j.err

# Bulk composition rebuild on the 774-unit corrected fragment counts.
# Estimator, reference and script are those behind accepted_composition.tsv.gz:
# Analysis/Deconvolution/scripts/12_run_bayesprism.R (InstaPrism, celltype labels) with
# Analysis/Deconvolution/reference/reference_human_scalesc_sce.rds, one run per cohort.
# A control run on the unchanged kallisto GSE126848 matrix must reproduce the stored proportions.
#   sbatch --export=ALL,COMP_ROOT=<new dir under Analysis/Deconvolution> run_composition_rebuild.sbatch
set -eo pipefail
umask 027

: "${COMP_ROOT:?COMP_ROOT is required}"
test ! -e "${COMP_ROOT}"
PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RNASEQ_ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
DECONV="${PROJECT_ROOT}/Analysis/Deconvolution"
REF_SCE="${DECONV}/reference/reference_human_scalesc_sce.rds"
SCRIPT12="${DECONV}/scripts/12_run_bayesprism.R"
test -f "${REF_SCE}" && test -f "${SCRIPT12}"

echo "step 1: corrected fragment counts to per-cohort symbol matrices"
"${RNASEQ_ENV}/bin/Rscript" --vanilla \
  "${PROJECT_ROOT}/scripts/analysis/bulk_refit_followup/prepare_corrected_composition_counts.R" \
  "${COMP_ROOT}/bulk"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate "${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv"
set -u
export R_LIBS_USER=""
R_PREFIX="${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv"
export R_LIBS_SITE="${R_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${R_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
R_BIN="${R_PREFIX}/bin/R"

run_prism() {  # counts_tsv result_dir dataset
  echo "step 2: InstaPrism $3 ($(date))"
  "${R_BIN}" --no-save --no-restore --slave --file="${SCRIPT12}" \
    --args "${REF_SCE}" "$1" "$2" "$3" celltype
}

for DATASET in GSE126848 GSE130970 GSE135251 GSE162694 GSE213621; do
  run_prism "${COMP_ROOT}/bulk/${DATASET}/${DATASET}_counts.tsv" "${COMP_ROOT}/results/${DATASET}" "${DATASET}"
done
# Reproduction control: the unchanged kallisto matrix behind the stored GSE126848 proportions.
run_prism "${DECONV}/bulk/GSE126848/GSE126848_counts.tsv" "${COMP_ROOT}/control_kallisto_input/GSE126848" GSE126848
echo "COMPOSITION_REBUILD_COMPLETE ${COMP_ROOT}"
