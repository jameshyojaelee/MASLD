#!/bin/bash
#SBATCH --job-name=1597
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1597_h3_recovery_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1597_h3_recovery_%j.err
set -euo pipefail
module load Bowtie2/2.5.4-linux-x86_64
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ $# -ne 1 ]]; then
    echo 'Require prospectively reviewed recovery processor SHA256' >&2
    exit 2
fi
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python codex_receiving_h3_reads_recovery.py --expected-source-sha256 "$1"
