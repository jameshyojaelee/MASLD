#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=90:00:00
# B3: PRESPEC_B tests 1-2. Mode from argument:
#   observed            -> one job
#   null (array 0-49)   -> 4 simulated datasets per task (200 total)
#   summary             -> decision table
# Usage: sbatch --job-name=<n> [--array=0-49] --output=<exec>/logs/%A_%a.out run_b3.sbatch <mode> <inputs.tsv.gz> <out_dir>
set -euo pipefail
mode=$1; inputs=$2; out=$3
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/aspect_latent
export PYTHONNOUSERSITE=1 XLA_FLAGS="--xla_force_host_platform_device_count=1" OMP_NUM_THREADS=4
mm=/gpfs/commons/home/jameslee/.local/bin/micromamba
if [ "$mode" = null ]; then
  "$mm" run -n allelic_qtl python b3_tests.py null --inputs "$inputs" --out "$out" --start $((SLURM_ARRAY_TASK_ID * 4)) --n 4
else
  "$mm" run -n allelic_qtl python b3_tests.py "$mode" --inputs "$inputs" --out "$out"
fi
