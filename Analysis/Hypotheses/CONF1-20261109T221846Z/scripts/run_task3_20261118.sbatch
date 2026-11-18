#!/bin/bash
#SBATCH --job-name=WS3A_T3
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/logs/ws3a_task3_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CONF1-20261109T221846Z/logs/ws3a_task3_%j.err
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
CONF1=$ROOT/Analysis/Hypotheses/CONF1-20261109T221846Z
XMOD2=$ROOT/Analysis/Hypotheses/XMOD2-20261113T151052Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

if [[ ! -f "$CONF1/out_c5_extract_20261118/pseudobulk_hep.tsv.gz" ]]; then
  echo "=== 1) c5_gse289173_extract_c1b_20261118.py ==="
  "$PY" -I "$CONF1/scripts/c5_gse289173_extract_c1b_20261118.py" "$CONF1/out_c5_extract_20261118"
fi
echo "=== 2) c5_trio_gse289173_20261118.R ==="
"$RS" "$CONF1/scripts/c5_trio_gse289173_20261118.R" "$CONF1/out_c5_extract_20261118" "$CONF1/out_c5_trio_20261118" TRIO
echo "=== 3) c2_methionine_hepatocyte_20261118.R ==="
"$RS" "$CONF1/scripts/c2_methionine_hepatocyte_20261118.R" "$CONF1/out_c2_20261118" TRIO CONF1_C2
echo "=== 4) XMOD2 c5_gse289173_xmod2_20261118.R ==="
"$RS" "$XMOD2/code/c5_gse289173_xmod2_20261118.R"
echo "=== 5) XMOD2 c2_hepatocyte_xmod2_20261118.R ==="
"$RS" "$XMOD2/code/c2_hepatocyte_xmod2_20261118.R"
echo "=== Task 3 Complete ==="
