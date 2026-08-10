#!/bin/bash
#SBATCH --job-name=studyaf
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=160G
#SBATCH --time=48:00:00
#SBATCH --array=1-35
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/study_af/studyaf_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/study_af/studyaf_%A_%a.err

# One study per task.  Each writes only its own fragment under
# config/af_sources/, so tasks cannot race on the shared override table; run
# `00b_build_study_af.R --merge` afterwards to assemble it.
#
# 160G because the join_lift studies (MVP x12, UKBB x3, FinnGen x2) hold a
# multi-million-row GRanges for liftOver alongside two large data.tables.
set -eo pipefail

PROJECT_ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
SRC="${PROJECT_ROOT}/GWAS/finemapping/src/seqfunc_v2/finemap"
FM="${PROJECT_ROOT}/GWAS/finemapping"

set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

# Study list = the studies of one placement, ordered exactly as the registry
# join orders them, so array index -> study is stable across resubmissions.
# AF_PLACEMENT=main (default, 35 studies, --array=1-35) is the finemap
# portfolio; AF_PLACEMENT=supp (15 studies, --array=1-15) is the tier-3/4
# remainder of the 50-study registry that canonical COLOC reads.  Selecting one
# placement at a time keeps the index stable and avoids re-running the
# expensive liftOver for studies already built.
AF_PLACEMENT="${AF_PLACEMENT:-main}"
STUDY=$(Rscript --vanilla -e '
suppressPackageStartupMessages(library(data.table))
fm <- commandArgs(TRUE)[1]; i <- as.integer(commandArgs(TRUE)[2]); pl <- commandArgs(TRUE)[3]
tiers <- fread(file.path(fm,"config","gwas_trait_tier.tsv"))
reg   <- fread(file.path(fm,"config","gwas_registry.tsv"))
m <- merge(tiers[placement==pl, .(study_name)], reg[, .(study_name)], by="study_name")
setorder(m, study_name)
cat(m$study_name[i])' "${FM}" "${SLURM_ARRAY_TASK_ID}" "${AF_PLACEMENT}")
[[ -n "${STUDY}" ]] || { echo "no study for placement=${AF_PLACEMENT} index ${SLURM_ARRAY_TASK_ID}" >&2; exit 2; }

echo "[$(date -u +%FT%TZ)] task ${SLURM_ARRAY_TASK_ID} -> ${STUDY}"
Rscript "${SRC}/00b_build_study_af.R" "${STUDY}"
echo "[$(date -u +%FT%TZ)] done ${STUDY}"
