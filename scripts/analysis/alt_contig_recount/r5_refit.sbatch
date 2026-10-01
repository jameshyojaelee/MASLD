#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
# R5: corrected-controls refit (the paper-candidate bulk model) on ALT-folded counts.
# Same environment and settings as scripts/manuscript/resource_f_five/run_refit_corrected_controls.sbatch.
# Usage: sbatch --job-name=<n> --output=<log> r5_refit.sbatch <corr_root> <output_root: .../bulk-corrected-controls-alt-UTC>
set -euo pipefail
umask 027
project_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
rnaseq_env=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export REFIT_ALT_ROOT=$1 REFIT_OUTPUT_ROOT=$2 REFIT_MODE=fit MASLD_PROJECT_ROOT=$project_root
test ! -e "$REFIT_OUTPUT_ROOT"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 PATH="$rnaseq_env/bin:$PATH"
cd "$project_root"
"$rnaseq_env/bin/Rscript" --vanilla "$project_root/scripts/analysis/alt_contig_recount/r5_refit_corrected_controls_alt.R"
echo "R5 complete $REFIT_OUTPUT_ROOT"
