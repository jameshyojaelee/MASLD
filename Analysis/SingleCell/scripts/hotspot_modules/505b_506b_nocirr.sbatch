#!/bin/bash
# CIRRHOSIS-EXCLUDED hepatocyte downstream (remediation R1), steps 2-3 of 3.
# Mirrors 505_phenotype_correlation.sbatch / 506_507.sbatch: cpu/32G/4cpu, rnaseq env.
# Runs phenotype correlation (505b) then bulk replication vs C2 canonical (506b),
# both writing into results_gpu_v2/hotspot_modules/hepatocytes_nocirrhosis/.
# Submitted with --dependency=afterok on the 501 nocirrhosis job.
#SBATCH --job-name=hotspot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/505b_506b_nocirr_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/505b_506b_nocirr_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

# Node-local micromamba proc lock (see note in 501_nocirrhosis_hep.sbatch) to
# avoid blocking on a cross-node ~/.cache/mamba/proc lock.
export XDG_CACHE_HOME="${SLURM_TMPDIR:-/tmp}/mambacache_${SLURM_JOB_ID}"
mkdir -p "$XDG_CACHE_HOME"
trap 'rm -rf "$XDG_CACHE_HOME"' EXIT

echo "[$(date)] 505b phenotype correlation (nocirrhosis)"
micromamba run -n rnaseq Rscript \
  Analysis/SingleCell/scripts/hotspot_modules/505b_phenotype_correlation_nocirr.R

echo "[$(date)] 506b bulk replication vs C2 canonical (nocirrhosis)"
micromamba run -n rnaseq Rscript \
  Analysis/SingleCell/scripts/hotspot_modules/506b_bulk_replication_nocirr.R

echo "[$(date)] Done 505b + 506b nocirrhosis"
