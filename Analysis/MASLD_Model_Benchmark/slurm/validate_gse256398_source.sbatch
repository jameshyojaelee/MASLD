#!/usr/bin/env bash
#SBATCH --job-name=source-validate-256398
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "${ROOT}"
bash -n slurm/acquire_audit_gse256398_h5.sbatch
module load h5py/3.9.0-foss-2023a
PYTHONPATH=${ROOT}:${PYTHONPATH:-} python -m unittest -v tests.unit.test_audit_gse256398_h5
sha256sum \
  config/datasets/gse256398.toml \
  scripts/audit_gse256398_h5.py \
  tests/unit/test_audit_gse256398_h5.py \
  slurm/acquire_audit_gse256398_h5.sbatch \
  slurm/validate_gse256398_source.sbatch
