#!/usr/bin/env bash
#SBATCH --job-name=model-check-311
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err
set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
ACTIVATION=${ROOT}/executions/pxd051911-activation-readiness-21109461
printf '%s  %s\n' "b2a5fea26985d7b10fe553afadbb693ea0bc47b6093d6567971cf4d2e5efade7" \
  "${ACTIVATION}/audit/readiness/participant_join.tsv" | sha256sum --check --strict
python scripts/compute_protein_transport_power_floors.py \
  --join "${ACTIVATION}/audit/readiness/participant_join.tsv" \
  --output "${ROOT}/config/campaigns/protein_transport_power_floors_20260825.json" \
  --draws 100000
