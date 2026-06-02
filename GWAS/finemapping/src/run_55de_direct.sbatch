#!/bin/bash
#SBATCH --job-name=gwasatac
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/run_55de_direct_%j.log
# ============================================================================
# run_55de_direct.sbatch -- finish the regulatory atlas (55d cCRE + 55e combine)
# using the rnaseq env's R binary DIRECTLY (PATH-prepended), bypassing the
# `micromamba run` environment-activation flock that is currently storm-
# contended by ~38 concurrent sessions. Idempotent with run_55bcd_regulatory.sbatch.
# ============================================================================
set -euo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC="$BASE/GWAS/finemapping/src"
export MASLD_PROJECT_ROOT="$BASE"
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export PATH="$ENV/bin:$PATH"
export LD_LIBRARY_PATH="$ENV/lib:${LD_LIBRARY_PATH:-}"
export R_HOME="$ENV/lib/R"
cd "$BASE"
echo "=== node: $(hostname) | start: $(date) | Rscript: $(which Rscript) ==="
echo "######## 55d SCREEN cCRE overlap (direct R) ########"
"$ENV/bin/Rscript" "$SRC/55d_ccre_overlap.R"
echo "######## 55e combine -> regulatory_atlas_columns.tsv (direct R) ########"
"$ENV/bin/Rscript" "$SRC/55e_combine_regulatory_atlas.R"
echo "=== done: $(date) ==="
