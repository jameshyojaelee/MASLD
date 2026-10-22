#!/bin/bash
#SBATCH --job-name=V6gates
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V6-20261021T234432Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V6-20261021T234432Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V6-20261021T234432Z/logs/%j.err
# V6 job 1: rebuild counts, gates K0, K3, K4. No program score.
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V6-20261021T234432Z
sha256sum "$RUN/PRESPEC_V6.md" "$RUN"/scripts/*.R > "$RUN/script_sha256_job1.txt"
(cd "$RUN/downloads_untrusted" && sha256sum -c SHA256SUMS) > "$RUN/logs/input_sha256_check_job1.txt"
/gpfs/commons/home/jameslee/micromamba/envs/geomx/bin/Rscript "$RUN/scripts/01_gates.R"
