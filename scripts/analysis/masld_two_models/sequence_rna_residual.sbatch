#!/bin/bash
#SBATCH --job-name=masld_seq_rna_residual
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_seq_rna_residual_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_seq_rna_residual_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=$PROJ/GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/sequence_rna_residual_${SLURM_JOB_ID}
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$PROJ/scripts/analysis/masld_two_models/sequence_rna_residual.py" --out "$OUT"
