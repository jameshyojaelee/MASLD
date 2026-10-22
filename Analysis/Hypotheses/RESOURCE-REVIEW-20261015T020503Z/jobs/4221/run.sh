#!/bin/bash
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/4221
mkdir -p "$J/tmp" "$J/fakeroot"
export TMPDIR="$J/tmp" XDG_CACHE_HOME="$J/tmp" MPLCONFIGDIR="$J/tmp"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
# Engine module creates an output dir under MASLD_PROJECT_ROOT at import; point it at a scratch root.
export MASLD_PROJECT_ROOT="$J/fakeroot"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd "$J"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python "$J/rebuild_observed.py"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -c "import numpy, pandas, anndata, scipy; print('numpy', numpy.__version__, 'pandas', pandas.__version__, 'anndata', anndata.__version__, 'scipy', scipy.__version__)"
