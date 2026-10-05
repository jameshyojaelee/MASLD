#!/bin/bash
#SBATCH --job-name=1584
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1584_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1584_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ $# -ne 2 ]]; then
    echo 'Require replacement primary job ID and reviewed primary source SHA256' >&2
    exit 2
fi
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python residual_pca_supplement.py --primary-job-id "$1" --expected-primary-source-sha256 "$2"
