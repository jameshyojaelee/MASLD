#!/bin/bash
#SBATCH --job-name=5801
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H58-20261015T112423Z/jobs/5801/slurm-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H58-20261015T112423Z/jobs/5801/slurm-%j.err
set -euo pipefail
R=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H58-20261015T112423Z
OUT=$R/jobs/5801/debug-$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}
mkdir "$OUT"; mkdir "$OUT/tmp"
export TMPDIR=$OUT/tmp XDG_CACHE_HOME=$OUT/tmp MPLCONFIGDIR=$OUT/tmp PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
cp $R/PRESPEC_H58.md "$OUT"/
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
/usr/bin/time -v $RS $R/scripts/02_reannotate_sites.R --out "$OUT"
/usr/bin/time -v $PY $R/scripts/01_gate_and_controls.py --out "$OUT"
$PY $R/scripts/03_estimates.py --out "$OUT" --debug-permute
