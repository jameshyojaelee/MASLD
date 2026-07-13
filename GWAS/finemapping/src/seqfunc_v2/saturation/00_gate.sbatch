#!/usr/bin/env bash
#SBATCH --job-name=saturation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/haplotype_saturation/v2/logs/gate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/haplotype_saturation/v2/logs/gate_%j.err
set -euo pipefail
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
OUT=${SATURATION_V2_ROOT:-$ROOT/GWAS/finemapping/results/seqfunc/haplotype_saturation/v2}
mkdir -p "$OUT/logs"
python3 "$ROOT/GWAS/finemapping/src/seqfunc_v2/saturation/00_upstream_gate.py" "$@"
