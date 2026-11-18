#!/bin/bash
#SBATCH --job-name=WS3A_T4
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/ws3a_task4_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/ws3a_task4_%j.err
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
DIG1=$ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

echo "=== Running 62b_genotype_table_20261118.py ==="
"$PY" -I "$DIG1/scripts/62b_genotype_table_20261118.py" "$DIG1/out_genotype_table_20261118"
echo "=== Running 63_genotype_state_20261118.R ==="
"$RS" "$DIG1/scripts/63_genotype_state_20261118.R"
echo "=== Task 4 Complete ==="
