#!/usr/bin/env bash
#SBATCH --job-name=masld-sealed-tamper
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1

module load python/3.11.5-GCCcore-13.2.0
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark

PYTHONPATH=src python -m unittest -v \
  tests.unit.test_tournament.SealedArtifactEnvelopeTamperTests.test_multiplicity_document_tamper_is_rejected_before_semantic_parse
