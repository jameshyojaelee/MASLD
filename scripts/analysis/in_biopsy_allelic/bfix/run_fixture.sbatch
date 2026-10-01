#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
# B-FIX synthetic fixture. Usage: sbatch --job-name=<n> --output=<out>/logs/%j.out run_fixture.sbatch <out>
set -euo pipefail
PY="env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
S=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/bfix/make_fixture.py
$PY $S "$1"
sha256sum $S > "$1/sha256_generator.txt"
(cd "$1" && find . -type f ! -path "./logs/*" ! -name "SHA256SUMS" -exec sha256sum {} + | sort -k2 > SHA256SUMS)
