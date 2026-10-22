#!/bin/bash
#SBATCH --job-name=V3score
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
# V3: lineage assignment and pseudobulk (Python), then frozen-weight scoring,
# HC3 models and random-set null (R).
# Usage: sbatch run_v3.sbatch <RESULT_DIR>
set -euo pipefail
RES="${1:?result dir required}"
ENVS=/gpfs/commons/home/jameslee/micromamba/envs
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/V3_program_holdout
export OMP_NUM_THREADS=8 PYTHONHASHSEED=0
cd "$RES"
(cd "$SCRIPTS" && sha256sum annotate_pseudobulk.py score_and_test.R) > "$RES/script_sha256.txt"
(cd downloads_untrusted && sha256sum -c SHA256SUMS)
"$ENVS/rapids_singlecell/bin/python" -I "$SCRIPTS/annotate_pseudobulk.py" "$RES"
"$ENVS/rapids_singlecell/bin/python" -I -m pip freeze > "$RES/pip_freeze.txt"
"$ENVS/rnaseq/bin/Rscript" "$SCRIPTS/score_and_test.R" "$RES"
