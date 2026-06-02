#!/bin/bash
#SBATCH --job-name=S4_hotspot_null
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/S4_hotspot_null_%j.out
#SBATCH --error=logs/S4_hotspot_null_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export HOTSPOT_NULL_PERM=${HOTSPOT_NULL_PERM:-100}
export HOTSPOT_NULL_CELLTYPE=${HOTSPOT_NULL_CELLTYPE:-hepatocytes}
export HOTSPOT_JOBS=${HOTSPOT_JOBS:-8}
export HOTSPOT_NEIGHBORS=${HOTSPOT_NEIGHBORS:-30}
# subsample for tractability: hepatocyte panel is ~657K cells
export HOTSPOT_NULL_SUBSAMPLE=${HOTSPOT_NULL_SUBSAMPLE:-60000}

# Ensure hotspot installed in rapids env (it's pure python, no GPU dep)
python -c "import hotspot" 2>/dev/null || pip install --quiet --user hotspotsc

python Analysis/SingleCell/scripts/hotspot_modules/502_empirical_null.py
