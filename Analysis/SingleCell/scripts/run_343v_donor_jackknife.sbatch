#!/bin/bash
#SBATCH --job-name=343v_jknf
#SBATCH --partition=gpu
#SBATCH --qos=interactive
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=Analysis/SingleCell/scripts/logs_scvi_loocv/343v_%x_%j.log
#SBATCH --error=Analysis/SingleCell/scripts/logs_scvi_loocv/343v_%x_%j.err
#
# Per-donor jackknife. Submit one job per donor with --export=DONOR=<sample_id>
# Example:
#   sbatch --export=ALL,DONOR=SRR19129274 --job-name=343v_F0a run_343v_donor_jackknife.sbatch

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation

source "$(micromamba info --base 2>/dev/null | awk -F': *' '/base environment/ {print $2}')/etc/profile.d/micromamba.sh" 2>/dev/null || true
eval "$(micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "[run] donor=${DONOR:?must export DONOR=<sample_id>}"
nvidia-smi || true

python Analysis/SingleCell/scripts/343v_scvi_donor_jackknife.py --hold-out-donor "${DONOR}"
