#!/bin/bash
#SBATCH --job-name=DIG1array
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z/logs/%j.err
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
sha256sum "$RUN/scripts/59_microarray_ballooning_replication.py" > "$RUN/script_sha256_job59.txt"
(cd "$RUN/downloads_untrusted/geo_microarray_components" && sha256sum -c SHA256SUMS)
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I "$RUN/scripts/59_microarray_ballooning_replication.py" "$RUN/downloads_untrusted/geo_microarray_components" "$RUN/out_microarray_ballooning"
