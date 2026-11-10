#!/bin/bash
#SBATCH --job-name=bf3_marker_split
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z/logs/%x_%j.out
set -euo pipefail
BF3=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/BF3-20261023T154841Z
BUILD=$BF3/build-20261023T155101Z
PRIMARY=$BF3/out-20261023T155821Z
OUT=$BF3/marker-split-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p $OUT
sha256sum $BF3/scripts/bf3_marker_split.py $BF3/DEVIATIONS.md > $OUT/script_sha256.txt
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
$PY -m pip freeze > $OUT/pip_freeze.txt
$PY $BF3/scripts/bf3_marker_split.py $BUILD $PRIMARY $OUT > $OUT/stdout.txt
echo done
