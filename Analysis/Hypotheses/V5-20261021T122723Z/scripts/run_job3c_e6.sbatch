#!/bin/bash
#SBATCH --job-name=V5e6
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z/logs/%j.err
# V5 job 3: pseudobulks (python), then scores, primary test, controls, descriptive arms (R).
# Both scripts stop unless CHECKED_BEFORE_OUTCOME.txt says PROCEED: yes.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=8
ENVS=/gpfs/commons/home/jameslee/micromamba/envs
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z
LIB=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/manuscript/program_observability_map/analysis_lib.R
sha256sum "$RUN/PRESPEC_V5.md" "$RUN/CHECKED_BEFORE_OUTCOME.txt" "$RUN"/scripts/*.py "$RUN"/scripts/*.R "$LIB" > "$RUN/script_sha256_job3c.txt"
"$ENVS/rnaseq/bin/Rscript" "$RUN/scripts/05_e6.R"
