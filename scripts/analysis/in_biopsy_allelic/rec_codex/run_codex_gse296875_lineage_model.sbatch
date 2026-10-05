#!/usr/bin/env bash
#SBATCH --job-name=1595
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1595_lineage_model_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1595_lineage_model_%j.err
set -euo pipefail
umask 027
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
exec /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$TASK_ROOT/scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse296875_lineage_model.py" \
  --expected-script-sha256 ccb6552d61b31f0d3fcf607f6b677c2b44909eaa6bdd6c1ae9777b2ec8e90430 \
  --execution-decision "$TASK_ROOT/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/20261001_codex_lineage_source_execution_decision.json" \
  --execution-decision-sha256 4fbe34ae9acf5b2c468d5f75fc59965810379f1a738504806bcb20a3c083e57f
