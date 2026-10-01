#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
# B5: PRESPEC_B test 3. Mode observed | null (array; 10 draws per task) | summary.
# Usage: sbatch --job-name=<n> [--array=0-99] --output=<exec>/logs/%A_%a.out run_b5.sbatch <mode> <inputs.tsv.gz> <b3_dir> <out_dir>
set -euo pipefail
mode=$1; inputs=$2; b3=$3; out=$4
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/aspect_latent
export PYTHONNOUSERSITE=1 XLA_FLAGS="--xla_force_host_platform_device_count=1" OMP_NUM_THREADS=4
mm=/gpfs/commons/home/jameslee/.local/bin/micromamba
if [ "$mode" = null ]; then
  "$mm" run -n allelic_qtl python b5_test3.py null --inputs "$inputs" --b3-dir "$b3" --out "$out" --start $((SLURM_ARRAY_TASK_ID * 10)) --n 10
else
  "$mm" run -n allelic_qtl python b5_test3.py "$mode" --inputs "$inputs" --b3-dir "$b3" --out "$out"
fi
