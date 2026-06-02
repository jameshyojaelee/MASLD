#!/bin/bash
# run_mouse_loo_dream.sbatch
# ---------------------------------------------------------------------------
# Refits the pooled mouse dream mega-analysis with one diet model held out at
# a time. Outputs feed the mouse |log2FC| stability analysis
# (M_loo_lfc_stability.R).
#
# Diet models are auto-discovered from config/mouse_datasets.yaml.
# The --array range should be set at submission time to match N diets:
#
# Submit (recommended — auto-computes array range):
#   cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#   N=$(python3 -c "import yaml; print(len(yaml.safe_load(open('config/mouse_datasets.yaml'))['diet_models'])-1)")
#   sbatch --array=0-${N}%4 RNA-seq/Mouse/Unified_Integration/scripts/run_mouse_loo_dream.sbatch
#
# Or with a fixed upper bound (excess tasks exit gracefully):
#   sbatch --array=0-19%4 RNA-seq/Mouse/Unified_Integration/scripts/run_mouse_loo_dream.sbatch
# ---------------------------------------------------------------------------
#SBATCH --job-name=mouse_loo_dream
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/Mouse/Unified_Integration/logs/loo_dream_%A_%a.log

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

# Auto-discover diet models from YAML config
YAML_FILE="config/mouse_datasets.yaml"
mapfile -t DIETS < <(python3 -c "
import yaml, sys
with open('$YAML_FILE') as f:
    cfg = yaml.safe_load(f)
diets = cfg.get('diet_models', [])
if not diets:
    print('ERROR: no diet_models in YAML', file=sys.stderr); sys.exit(1)
for d in diets:
    print(d)
")

N_DIETS=${#DIETS[@]}
if [[ $N_DIETS -eq 0 ]]; then
  echo "ERROR: no diet models found in $YAML_FILE" >&2
  exit 1
fi

# Gracefully exit if this array task index exceeds the number of diets
if [[ $SLURM_ARRAY_TASK_ID -ge $N_DIETS ]]; then
  echo "[$(date)] Task $SLURM_ARRAY_TASK_ID >= $N_DIETS diets; nothing to do. Exiting."
  exit 0
fi

# Map array index -> diet model
export LOO_DIET="${DIETS[$SLURM_ARRAY_TASK_ID]}"

echo "[$(date)] Diet models ($N_DIETS): ${DIETS[*]}"
echo "[$(date)] Running LOO fold: $LOO_DIET (task $SLURM_ARRAY_TASK_ID of 0-$((N_DIETS-1)))"

Rscript RNA-seq/Mouse/Unified_Integration/scripts/M_loo_dream_refits.R

echo "[$(date)] Done: $LOO_DIET"
