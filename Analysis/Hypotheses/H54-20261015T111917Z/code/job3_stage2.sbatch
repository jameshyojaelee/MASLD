#!/bin/bash
# Job 3: stage 2 (outcome). Submitted by hand only after stage1/stop_rule.tsv is read.
#SBATCH --job-name=5403
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/logs/%x_%j.err
# env: /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript; OMP/OPENBLAS threads 4
# Real run: no H54_* variables set; writes <run>/stage2 and <run>/RESULT_H54.md.
# Debug (gene-permuted only): H54_S2_OUT, H54_S1_DIR, H54_DEBUG_GENE_PERMUTED=1, H54_DEBUG_DIVISOR, H54_DEBUG_FORCE_OPEN_SECONDARY.
set -euo pipefail
umask 027
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z
OUT="${H54_S2_OUT:-${RUN}/stage2}"
S1="${H54_S1_DIR:-${RUN}/stage1}"
test ! -e "${OUT}"
test -s "${S1}/stop_rule.tsv"
export TMPDIR="${RUN}/tmp_${SLURM_JOB_ID}"; mkdir "${TMPDIR}"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 FLEXIBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla "${RUN}/code/stage2.R" "${OUT}" "${S1}"
echo "H54_STAGE2_COMPLETE ${OUT} ($(date))"
