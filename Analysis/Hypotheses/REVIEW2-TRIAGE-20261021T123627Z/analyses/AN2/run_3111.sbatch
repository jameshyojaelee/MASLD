#!/bin/bash
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN2/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN2/%j.out
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
N=$ROOT/figures/candidates/fig2d-gws-leads-20261021T124551Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR=$N/logs MPLCONFIGDIR=$N/logs
cd "$ROOT"
$PY -m pip freeze > $N/logs/pip_freeze.txt
sha256sum $N/code/*.py $N/code/*.R > $N/logs/code_sha256.txt
$PY $N/code/compute.py --out $N/values
for mode in gws_p gws_tier; do
  P=$N/panels_$mode
  mkdir -p $P/code/originals $P/logs $P/check
  cp -p $N/code/originals/*.R $P/code/originals/
  $PY $N/code/build_tables.py $N/values/$mode $P
  $RS $N/code/fig2D_composition.R $P
  $RS $N/code/fig2D_ancestry.R $P
  for f in Fig2D_coding_mass_by_trait_class Fig2D_coding_by_ancestry; do
    pdftotext -layout $P/$f.pdf $P/check/$f.txt
    pdftoppm -png -r 300 -singlefile $P/$f.pdf $P/check/$f
  done
done
echo "AN2 DONE $(date -u)"
