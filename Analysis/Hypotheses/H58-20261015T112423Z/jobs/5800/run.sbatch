#!/bin/bash
#SBATCH --job-name=5800
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800/slurm-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800/slurm-%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
JOBDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800
OUT=$JOBDIR/H58-$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}
mkdir "$OUT"                      # fails if it exists; nothing is overwritten
mkdir -p "$OUT/tmp" "$OUT/scripts"
export TMPDIR=$OUT/tmp XDG_CACHE_HOME=$OUT/tmp MPLCONFIGDIR=$OUT/tmp PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
cp "$JOBDIR"/PRESPEC_H58.md "$JOBDIR"/inputs.sha256 "$JOBDIR"/inputs_g7.sha256 "$OUT"/
cp "$JOBDIR"/scripts/* "$OUT/scripts/"
sha256sum "$OUT"/scripts/* > "$OUT/script_sha256.txt"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
$PY -m pip freeze > "$OUT/pip_freeze.txt"
$RS -e 'sessionInfo()' > "$OUT/sessionInfo.txt"
cd "$ROOT"
if ! sha256sum -c "$OUT/inputs.sha256" "$OUT/inputs_g7.sha256" > "$OUT/input_check.txt" 2>&1; then
  printf 'G0 failed: input sha256 mismatch (see input_check.txt). Nothing computed.\n' > "$OUT/RESULT_H58.md"
  exit 1
fi
echo "G0 PASS: all input sha256 match"
$RS "$OUT/scripts/02_reannotate_sites.R" --out "$OUT"
$PY "$OUT/scripts/01_gate_and_controls.py" --out "$OUT"
$PY "$OUT/scripts/03_estimates.py" --out "$OUT"
$PY "$OUT/scripts/04_fig2a_support.py" --out "$OUT"
echo "H58 DONE"
