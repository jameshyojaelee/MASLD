#!/bin/bash
#SBATCH --job-name=343w_array
#SBATCH --partition=gpu
#SBATCH --gres=gpu:b6k:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --qos=nslab
#SBATCH --array=0-6
#SBATCH --output=Analysis/SingleCell/scripts/logs_scvi_loocv/343w_%a_%j.log
#SBATCH --error=Analysis/SingleCell/scripts/logs_scvi_loocv/343w_%a_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p Analysis/SingleCell/scripts/logs_scvi_loocv
eval "$(micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell
set -u

COHORTS=(GSE136103 GSE174748 GSE185477 GSE189600 GSE202379 GSE244832 Liver_Atlas)
COHORT=${COHORTS[$SLURM_ARRAY_TASK_ID]}
echo "[343w] holding out: $COHORT"
python Analysis/SingleCell/scripts/343w_scvi_cohort_out.py --hold-out-cohort "$COHORT"
