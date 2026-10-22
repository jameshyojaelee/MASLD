#!/bin/bash
#SBATCH --job-name=h50_plasma_panel
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H50-20261009T135841Z/logs/slurm-%x-%j.out
# H50: PXD051911 count recheck, detection gate, positive control. Single job, no array.
# Writes Analysis/Hypotheses/H50-20261009T135841Z/results/run-<UTC>/ (new directory).
set -euo pipefail
hyp=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H50-20261009T135841Z
scr=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H50_plasma_panel
rscript=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run="${hyp}/results/run-$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "${run}"
echo "host $(hostname) job ${SLURM_JOB_ID} run ${run}"
python3 --version > "${run}/python_version.txt"
python3 -I -m pip freeze > "${run}/pip_freeze.txt" 2>&1 || true
python3 -I "${scr}/01_counts_and_gate.py" "${run}/counts_and_gate"
"${rscript}" "${scr}/02_positive_control.R" "${run}/positive_control"
sha256sum "${scr}"/01_counts_and_gate.py "${scr}"/02_positive_control.R "${hyp}/PRESPEC_H50.md" > "${run}/script_and_prespec.sha256"
echo done
