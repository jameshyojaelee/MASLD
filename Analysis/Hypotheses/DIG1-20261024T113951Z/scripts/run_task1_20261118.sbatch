#!/bin/bash
#SBATCH --job-name=WS3A_T1
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/ws3a_task1_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/ws3a_task1_%j.err
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
DIG1=$ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z
CCL1=$ROOT/Analysis/Hypotheses/CCL1-20261110T235204Z

echo "=== Running 56_atlas_state_replication_20261118.py ==="
"$PY" -I "$DIG1/scripts/56_atlas_state_replication_20261118.py" "$DIG1/out_atlas_state_replication_20261118"
echo "=== Running 57_state_bulk_source_20261118.py ==="
"$PY" -I "$DIG1/scripts/57_state_bulk_source_20261118.py" "$DIG1/out_state_bulk_source_20261118"
echo "=== Running ccl1_identity_20261118.py ==="
"$PY" -I "$CCL1/scripts/ccl1_identity_20261118.py" "$CCL1/out_ccl1_20261118"
echo "=== Task 1 Complete ==="
