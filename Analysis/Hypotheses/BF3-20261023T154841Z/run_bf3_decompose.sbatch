#!/bin/bash
#SBATCH --job-name=bf3_decompose
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z/logs/%x_%j.out
set -euo pipefail
BF3=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z
BUILD=$BF3/build-20261023T155101Z
OUT=$BF3/out-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p $OUT
cd $BF3 && sha256sum -c PRESPEC_BF3.md.sha256
sha256sum $BF3/scripts/bf3_check_build.py $BF3/scripts/bf3_decompose.py > $OUT/script_sha256.txt
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY -m pip freeze > $OUT/pip_freeze.txt
$PY $BF3/scripts/bf3_check_build.py $BUILD > $OUT/check_build.txt
cat $OUT/check_build.txt
$PY $BF3/scripts/bf3_decompose.py $BUILD $OUT > $OUT/results_stdout.txt
echo done
