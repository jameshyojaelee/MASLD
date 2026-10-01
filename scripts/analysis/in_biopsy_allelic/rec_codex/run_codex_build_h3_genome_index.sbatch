#!/bin/bash
#SBATCH --job-name=1594
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
set -euo pipefail
module load Bowtie2/2.5.4-linux-x86_64
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex
export CODEX_REC_OUTPUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python codex_build_h3_genome_index.py --output "$CODEX_REC_OUTPUT/codex_h3_genome_index_${SLURM_JOB_ID}"
