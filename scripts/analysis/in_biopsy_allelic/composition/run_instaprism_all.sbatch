#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=24:00:00
# Usage: sbatch --job-name=<n> --output=<out>/%j.out run_instaprism_all.sbatch <out_dir>
set -euo pipefail
export COMP_OUT=$1
MD=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/music_deconv
R_LIBS_USER="" R_LIBS_SITE="$MD/lib/R/library" LD_LIBRARY_PATH="$MD/lib:${LD_LIBRARY_PATH:-}" \
  "$MD/bin/R" --no-save --no-restore --slave --file=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/composition/instaprism_all.R
