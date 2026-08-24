#!/usr/bin/env bash
#SBATCH --job-name=masld-sealed-promotion
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1

module load python/3.11.5-GCCcore-13.2.0
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark

PYTHONPATH=src python -m unittest -v \
  tests.unit.test_tournament.FrozenTournamentChainTests.test_rederived_metrics_holm_family_and_positive_promotion \
  tests.unit.test_tournament.FrozenTournamentChainTests.test_opened_outcome_is_terminal_but_never_promoted_without_metric_bundle \
  tests.unit.test_tournament.FrozenTournamentChainTests.test_sealed_metric_and_multiplicity_tamper_fail_closed
