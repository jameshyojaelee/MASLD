#!/bin/bash
#SBATCH --job-name=agp_ad1_build
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad1_build_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad1_build_%j.err

set -uo pipefail
export PYTHONNOUSERSITE=1
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba

MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
ENV=alphagenome_local
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=$PROJ/scripts/analysis/alphagenome_program
RUN=$PROJ/GWAS/finemapping/results/alphagenome_program/ad1-gate-20260915T114615Z

mkdir -p "$RUN/inputs" "$RUN/tables" "$RUN/logs" "$RUN/env" "$RUN/raw" "$RUN/sbatch"

echo "=== node $(hostname) job $SLURM_JOB_ID $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
$MM run -n $ENV python "$SRC/i1_00_build_gate_set.py" --out "$RUN"
echo "  exit=$?"
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
