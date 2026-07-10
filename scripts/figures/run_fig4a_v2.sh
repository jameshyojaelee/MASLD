#!/bin/bash
# Render the 5 Fig 4A validation-overview v2 candidates (fig4a_v2_render_all.py).
# job-name = matplotlib (per project SLURM naming: single tool word). rnaseq env
# (pandas/numpy; pure-numpy null, no scipy). matplotlib import can hang on the login
# node -> always render via SLURM. --output goes to SHARED gpfs (NOT /scratch, which
# is node-local/invisible from the login node); MPLCONFIGDIR can be node-local /scratch.
#SBATCH --job-name=matplotlib
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4a_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4a_v2_%j.err
# NOTE: do NOT `set -u` here — the rnaseq env's binutils activation script
# references ADDR2LINE unbound and aborts under set -u (known rnaseq quirk).
set -e

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u  # safe to enable strict mode AFTER the conda activation scripts have run

export MPLCONFIGDIR=/scratch/claude-91825/-gpfs-commons-groups-sanjana-lab-Cas13-MASLD-library-design/b8c2f7b0-a0b6-40c6-b7e9-2514eeb77901/scratchpad/mplcache
mkdir -p "$MPLCONFIGDIR"

# Pin hash randomization so the permutation null is reproducible across runs.
# permutation_null() iterates over Python sets when building its sampling index
# arrays; set-iteration order (hence the seed-42 rng.choice draws) otherwise
# varies per process, wobbling exp/p slightly (e.g. 51.4<->51.5, p 1e-3<->2e-3).
# PYTHONHASHSEED=0 makes string hashing deterministic -> stable null -> stable numbers.
export PYTHONHASHSEED=0

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures
python fig4a_v2_render_all.py
