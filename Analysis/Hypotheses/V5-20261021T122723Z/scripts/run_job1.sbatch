#!/bin/bash
#SBATCH --job-name=V5qc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z/logs/%j.err
# V5 job 1: G0-G3, labels, PC1, PC2, diagnostics. No program score.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
sha256sum "$RUN/PRESPEC_V5.md" "$RUN"/scripts/*.py > "$RUN/script_sha256_job1.txt"
(cd "$RUN/downloads_untrusted" && sha256sum -c SHA256SUMS) > "$RUN/logs/input_sha256_check_job1.txt"
"$PY" -I "$RUN/scripts/01_qc_annotate.py"
"$PY" -I -m pip freeze > "$RUN/pip_freeze_spatial.txt"
