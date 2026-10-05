#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
# C1b acquisition check for the Model C external test (Amendment 6): TCGA-LIHC non-tumor
# slides in IDC matched to RNA-seq sample barcodes; no image download, no RNA value read.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out run_c1b.sbatch <exec_dir>
set -euo pipefail
exec_dir=$1
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mm=/gpfs/commons/home/jameslee/.local/bin/micromamba
export PYTHONNOUSERSITE=1
mkdir -p "$exec_dir"
"$mm" run -n he_molecular python -m pip freeze > "$exec_dir/pip_freeze.txt"
"$mm" run -n he_molecular python "$repo/scripts/analysis/he_molecular/c1b_tcga_acquisition_check.py" --out "$exec_dir"
