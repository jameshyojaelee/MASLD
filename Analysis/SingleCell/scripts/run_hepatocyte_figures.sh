#!/bin/bash
# SBATCH wrapper for Script 311 (figures).
# Run manually AFTER reviewing subtype labels.
#
# Usage: sbatch run_hepatocyte_figures.sh

#SBATCH --job-name=hep_311_fig
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/hepatocyte_subclustering/311_fig_%j.out
#SBATCH --error=logs/hepatocyte_subclustering/311_fig_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
micromamba run -n rnaseq Rscript "${SCRIPT_DIR}/311_hepatocyte_subcluster_figures.R"
