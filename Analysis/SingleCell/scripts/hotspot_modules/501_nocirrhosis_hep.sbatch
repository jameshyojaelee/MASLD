#!/bin/bash
# CIRRHOSIS-EXCLUDED hepatocyte Hotspot run (remediation R1), step 1 of 3.
# Mirrors 501_run_hotspot.sbatch (the canonical hepatocyte run): bigmem/512G/16cpu,
# nslab QOS, rapids_singlecell env, HOTSPOT_MAX_GENES=2000. The 657k-cell
# hepatocyte subset needs ~512G, which exceeds the cpu partition's 210G cap, so
# this heavy step uses bigmem (exactly the proven canonical config). plain nslab
# QOS is used (NOT --qos=interactive) because the 4 interactive slots are taken
# and the nslab 7TB cap has ~2.3TB headroom. The light downstream R steps
# (505b, 506b) run on cpu in the chained companion script.
#
# --exclude-stage Cirrhosis filters at the DONOR/STAGE level (GSE202379 carries
# both cirrhosis and non-cirrhosis donors) and auto-routes output to
# results_gpu_v2/hotspot_modules/hepatocytes_nocirrhosis/ (never touches the
# canonical hepatocytes/ dir).
#SBATCH --job-name=hotspot
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/501_nocirr_hep_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/501_nocirr_hep_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs

# Relocate the micromamba proc lock to a node-local tmp dir. The default lock at
# ~/.cache/mamba/proc is shared across nodes on this filesystem; a long-running
# login-node micromamba process (e.g. a streamlit app) can otherwise make every
# `micromamba run` here block indefinitely on "Waiting for other mamba process".
export XDG_CACHE_HOME="${SLURM_TMPDIR:-/tmp}/mambacache_${SLURM_JOB_ID}"
mkdir -p "$XDG_CACHE_HOME"
trap 'rm -rf "$XDG_CACHE_HOME"' EXIT

export HOTSPOT_MAX_GENES=2000
export HOTSPOT_JOBS=16

echo "[$(date)] 501 hepatocytes --exclude-stage Cirrhosis (-> hepatocytes_nocirrhosis/)"
micromamba run -n rapids_singlecell python \
  Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py \
  --cell-type hepatocytes \
  --exclude-stage Cirrhosis

echo "[$(date)] Done 501 nocirrhosis hepatocytes"
