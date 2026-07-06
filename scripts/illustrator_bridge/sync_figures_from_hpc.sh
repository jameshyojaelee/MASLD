#!/usr/bin/env bash
# sync_figures_from_hpc.sh
# Pull the latest publication figure PDFs from the MASLD HPC project down to a
# local mirror on this Mac.
#
# RUN ON THE MAC (not the HPC). Requires: rsync (ships with macOS) + ssh access
# to the HPC login node. The Mac PULLS from the HPC (Macs are usually behind NAT,
# so the HPC can't push to you).
#
# Usage:
#   ./sync_figures_from_hpc.sh
# Override any variable inline, e.g.:
#   HPC_HOST=masld-hpc LOCAL_DIR=~/masld_figures/ ./sync_figures_from_hpc.sh
set -euo pipefail

# ---- CONFIG (edit these, or export before running) --------------------------
# SSH target for the HPC LOGIN node (use whatever you `ssh` into normally).
# Recommended: add a Host alias "masld-hpc" to ~/.ssh/config, then leave this.
HPC_HOST="${HPC_HOST:-masld-hpc}"

# Absolute path to figures/ on the HPC (do not change unless the repo moves).
HPC_FIG_DIR="${HPC_FIG_DIR:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/}"

# Local mirror on this Mac. Keep this SEPARATE from your .ai master documents.
LOCAL_DIR="${LOCAL_DIR:-$HOME/masld_figures/}"
# -----------------------------------------------------------------------------

mkdir -p "$LOCAL_DIR"

# Pull PDFs + the manifest + READMEs only (keeps the mirror small).
# --delete keeps the mirror faithful to the HPC. Excluded local files (e.g. any
# .ai you accidentally drop here) are PROTECTED from deletion (no --delete-excluded).
rsync -avz --delete --prune-empty-dirs \
  --include='*/' \
  --include='*.pdf' \
  --include='figures_manifest.tsv' \
  --include='README.md' \
  --exclude='*' \
  "${HPC_HOST}:${HPC_FIG_DIR}" "$LOCAL_DIR"

echo "Synced HPC figures -> $LOCAL_DIR"
