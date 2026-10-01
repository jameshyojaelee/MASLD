#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
# B2: Model B development inputs (program scores + labels), corrected and original counts.
# Usage: sbatch --job-name=<n> --dependency=afterok:<recounts> --output=<exec>/logs/%j.out run_b2.sbatch <multi_counts_dir> <exec_dir>
set -euo pipefail
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$repo/scripts/analysis/aspect_latent"
sha256sum -c PRESPEC.sha256
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript b2_build_inputs.R "$1" "$2"
