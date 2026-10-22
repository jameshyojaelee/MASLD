#!/bin/bash
#SBATCH --job-name=V0cov4
#SBATCH --constraint=v4
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V0-20261021T122533Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V0-20261021T122533Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V0-20261021T122533Z/logs/%j.err
# V0: rerun V3 labelling/pseudobulk, check reproduction, compute program-gene
# detection under the H90 rule. No outcome is computed.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=8
ENVS=/gpfs/commons/home/jameslee/micromamba/envs
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V0-20261021T122533Z
V3S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/V3_program_holdout
sha256sum "$RUN/PRESPEC_V0.md" "$RUN/scripts/v0_coverage_audit.py" "$V3S/annotate_pseudobulk.py" > "$RUN/script_sha256_v4node.txt"
(cd "$RUN/downloads_untrusted" && sha256sum -c SHA256SUMS) > "$RUN/logs/input_sha256_check_v4node.txt"
"$ENVS/rapids_singlecell/bin/python" -I "$RUN/scripts/v0_coverage_audit.py" "$RUN" out_v4node
"$ENVS/rapids_singlecell/bin/python" -I -m pip freeze > "$RUN/pip_freeze_v4node.txt"
