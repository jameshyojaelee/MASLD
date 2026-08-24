#!/usr/bin/env bash
#SBATCH --job-name=masld-gse296875-view
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1

module load python/3.11.5-GCCcore-13.2.0
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark

PYTHONPATH=src python -m masld_bench.cli registry validate --config-root config
PYTHONPATH=src python - <<'PY'
from masld_bench.registry import Registry

registry = Registry.load("config")
view = registry.get_dataset_view("gse296875_rna_atac_smoke_1000_v1")
assert view.parent_dataset_id == "gse296875"
assert view.row_count == 1000
assert view.biological_unit_count == 39
assert tuple(view.modalities) == ("single_nucleus_rna", "single_nucleus_atac")
assert tuple(item.value for item in view.pairing_levels) == ("same_nucleus",)
assert sum(view.class_counts.values()) == 1000
print(view.to_dict())
PY
