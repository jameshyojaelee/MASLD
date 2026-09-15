#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-check-040
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-040-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/model-cpu-check-040-%j.err

set -euo pipefail
umask 027

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
DRIVER=${ROOT}/slurm/test_and_freeze_scbasset_crossfold_evaluator.sbatch
printf '%s  %s\n' c15ba471c16b849f3fc348c2156e49b3cb7ea67d179faefda88d4a7fbb17a11d "${DRIVER}" \
  | sha256sum --check --strict

python() {
  if [[ "${1:-}" == -m && "${2:-}" == unittest ]]; then
    /gpfs/commons/home/jameslee/.local/bin/micromamba run \
      --prefix /gpfs/commons/home/jameslee/micromamba/envs/immune-atlas \
      env PYTHONNOUSERSITE=1 "PYTHONPATH=${PYTHONPATH:-}" python "$@"
  else
    command python "$@"
  fi
}
export -f python
source "${DRIVER}"
