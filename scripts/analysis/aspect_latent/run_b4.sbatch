#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# B4: final NUTS fit at the kept K (4 chains in parallel, one per CPU).
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out run_b4.sbatch <inputs.tsv.gz> <out_dir> [extra args]
set -euo pipefail
inputs=$1; out=$2; shift 2
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/aspect_latent
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1
/gpfs/commons/home/jameslee/.local/bin/micromamba run -n allelic_qtl python b4_nuts.py --inputs "$inputs" --out "$out" "$@"
