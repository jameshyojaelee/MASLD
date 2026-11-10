#!/bin/bash
#SBATCH --job-name=A2_stats
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage/logs/stats_%j.out
# Usage: sbatch run_04_09.sbatch [test]   (test: no LD scores, results to results_test/)
set -euo pipefail
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 TMPDIR=$A/work
E=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
if [[ "${1:-}" == "test" ]]; then
  export A2_ALLOW_NO_LD=1 A2_RES=$A/results_test
  mkdir -p $A/results_test && cp $A/results/q2_expression_lineage.tsv $A/results_test/
fi
for s in 04_annotate 05_q1 07_q2 08_q3 09_q4; do
  echo "== $s $(date -u +%H:%M:%S)"; $E/python $A/scripts/$s.py
done
if [[ "${1:-}" != "test" ]]; then
  echo "== sha256 $(date -u +%H:%M:%S)"; bash $A/scripts/10_sha256.sh
  $E/python -m pip freeze > $A/pip_freeze.txt
fi
echo "== done $(date -u +%H:%M:%S)"
