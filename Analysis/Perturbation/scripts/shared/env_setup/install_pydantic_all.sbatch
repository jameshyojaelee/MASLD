#!/bin/bash
#SBATCH --job-name=pydantic-inst
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/install_pydantic_%j.log

# Install pydantic into the envs that need it for runner_base output_schema.
# Use conda (not pip) to avoid user-site leak issue.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

ENVS_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs
eval "$(micromamba shell hook --shell bash)"

for E in perturbation_gears perturbation_tahoe perturbation_scgen_cpa perturbation_ccc; do
  ENV_PATH="$ENVS_DIR/$E"
  if [ ! -d "$ENV_PATH" ]; then
    echo "[skip] $E missing"; continue
  fi
  echo "==================== $E ===================="
  echo "[+] installing pydantic via conda"
  micromamba install -p "$ENV_PATH" -c conda-forge -y pydantic >/dev/null 2>&1 && \
    echo "  conda install OK" || echo "  conda install FAILED"
  # Verify
  micromamba run -p "$ENV_PATH" env PYTHONNOUSERSITE=1 python -c "
import pydantic, sys
print(f'  pydantic {pydantic.__version__} from {pydantic.__file__[:80]}')
print(f'  sys.path[0]: {sys.path[0]}')
"
done
echo "[done] pydantic install pass complete"
