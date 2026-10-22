#!/bin/bash
#SBATCH --job-name=NNNN
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H56-20261015T112556Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H56-20261015T112556Z/logs/%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RSCRIPT=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export H56_RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H56-20261015T112556Z
export H56_MODE=${H56_MODE:-final}
if [ "$H56_MODE" = final ]; then export H56_OUT=$H56_RUN/out; else export H56_OUT=$H56_RUN/debug/$SLURM_JOB_ID; mkdir "$H56_OUT"; fi
export TMPDIR=$H56_RUN/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 XDG_CACHE_HOME=$H56_RUN/tmp/cache
echo "start $(date -u) host $(hostname) mode $H56_MODE out $H56_OUT"
set +e
$PY $H56_RUN/scripts/h56_prep.py; rc1=$?
echo "prep exit $rc1"
rc2=NA
if [ $rc1 -eq 0 ]; then $RSCRIPT $H56_RUN/scripts/h56_fit.R; rc2=$?; echo "fit exit $rc2"; fi
$RSCRIPT $H56_RUN/scripts/h56_result.R
if [ -f $H56_OUT/per_gene_all_cells.tsv.gz ]; then $RSCRIPT $H56_RUN/scripts/h56_figure.R; echo "figure exit $?"; fi
$PY -m pip freeze > $H56_OUT/pip_freeze.txt 2>/dev/null
echo "end $(date -u)"
