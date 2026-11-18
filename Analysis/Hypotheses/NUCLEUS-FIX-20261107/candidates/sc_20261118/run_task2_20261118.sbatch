#!/bin/bash
#SBATCH --job-name=WS3A_T2
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/NUCLEUS-FIX-20261107/candidates/sc_20261118/ws3a_task2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/NUCLEUS-FIX-20261107/candidates/sc_20261118/ws3a_task2_%j.err
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
SC2=$ROOT/Analysis/Hypotheses/NUCLEUS-FIX-20261107/candidates/sc_20261118
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

if [[ ! -f "$SC2/fig4c_program_match_counts_bootstrap.tsv" ]]; then
  echo "=== Running run_dedup_lineage_specificity_fig4c.py ==="
  "$PY" -I "$SC2/run_dedup_lineage_specificity_fig4c.py"
fi
echo "=== Running run_dedup_pseudobulk_de_and_stouffer.R ==="
"$RS" "$SC2/run_dedup_pseudobulk_de_and_stouffer.R"
echo "=== Task 2 Complete ==="
