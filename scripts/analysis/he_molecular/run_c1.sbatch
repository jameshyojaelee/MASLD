#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
# C1 acquisition check for Model C: GTEx liver slides in the NCI Imaging Data
# Commons matched to v8 liver RNA samples; no image download, no RNA value read.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out run_c1.sbatch <exec_dir>
set -euo pipefail
exec_dir=$1
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mm=/gpfs/commons/home/jameslee/.local/bin/micromamba
export PYTHONNOUSERSITE=1
mkdir -p "$exec_dir"
"$mm" run -n he_molecular python -m pip freeze > "$exec_dir/pip_freeze.txt"
"$mm" run -n he_molecular python "$repo/scripts/analysis/he_molecular/c1_acquisition_check.py" \
  --sample-attributes "$repo/data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt" \
  --out "$exec_dir"
