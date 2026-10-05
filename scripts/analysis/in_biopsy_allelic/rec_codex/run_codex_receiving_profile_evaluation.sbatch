#!/bin/bash
#SBATCH --job-name=1598
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1598_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1598_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ $# -ne 2 ]]; then
    echo 'Require exact reviewed H3 recovery job ID and source SHA256' >&2
    exit 2
fi
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python codex_receiving_profile_evaluation.py --h3-recovery-job-id "$1" --expected-recovery-source-sha256 "$2"
