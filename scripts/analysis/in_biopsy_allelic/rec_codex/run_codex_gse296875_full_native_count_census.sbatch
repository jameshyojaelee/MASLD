#!/usr/bin/env bash
#SBATCH --job-name=1594
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1594_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1594_%j.err
set -euo pipefail
umask 027
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load R/4.4.3
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
exec Rscript --vanilla "$TASK_ROOT/scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse296875_full_native_count_census.R"
