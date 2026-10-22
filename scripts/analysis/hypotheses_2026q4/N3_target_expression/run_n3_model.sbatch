#!/bin/bash
#SBATCH --job-name=N3_model
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N3-20261009T140326Z/logs/%x-%j.out
# N3 step 3: donor-level limma-voom models, controls, permutations, bulk replication.
# Usage: sbatch run_n3_model.sbatch <OUT_SUBDIR> <N_PERM>
set -euo pipefail
OUTSUB="${1:?out subdir}"; NPERM="${2:?n perm}"
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RES=$REPO/Analysis/Hypotheses/N3-20261009T140326Z
SCR=$REPO/scripts/analysis/hypotheses_2026q4/N3_target_expression
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
test -f "$RES/donor_meta.tsv"
sha256sum "$SCR/03_donor_model.R" > "$RES/script_sha256_model_${OUTSUB}.txt"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$SCR/03_donor_model.R" "$RES" "$OUTSUB" "$NPERM" "${SLURM_CPUS_PER_TASK:-8}"
