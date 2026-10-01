#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# B1 design-only checks for Model B (NAS-sum unit test, K = 2 recovery, GSE281797
# minimum detectable effects). Reads no RNA. Job name and log set on the command line.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out run_b1.sbatch <exec_dir>
set -euo pipefail
exec_dir=$1
src=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/aspect_latent
mm=/gpfs/commons/home/jameslee/.local/bin/micromamba
export PYTHONNOUSERSITE=1
mkdir -p "$exec_dir"
"$mm" run -n allelic_qtl python -m pip freeze > "$exec_dir/pip_freeze.txt"
cd "$src"
"$mm" run -n allelic_qtl python b1_simulate.py --out "$exec_dir" --chains 4 --warmup 500 --draws 500 --sims 2000
