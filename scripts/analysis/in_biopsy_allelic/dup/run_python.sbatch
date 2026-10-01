#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
# Job name and log path are set on the sbatch command line (numbers only).
#
# Runs one B-Dup Python step (01_regions.py or 04_collect.py) in the rnaseq env and
# records the environment (pip freeze, samtools version) in <versions_dir> if not
# already there.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out run_python.sbatch <versions_dir> <script.py> [args...]
set -euo pipefail
versions=$1; shift
env=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
export PYTHONNOUSERSITE=1
[ -e "$versions/pip_freeze.txt" ] || "$env/python" -m pip freeze > "$versions/pip_freeze.txt"
[ -e "$versions/samtools_version.txt" ] || "$env/samtools" --version | head -2 > "$versions/samtools_version.txt"
"$env/python" "$@"
