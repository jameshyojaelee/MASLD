#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/coloc23_chain_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/coloc23_chain_%j.out
# Phase 4 of the 23-GWAS COLOC refactor — full downstream chain on healthy io node.
# Order: 27a (atlas, clean coloc layers) -> 77b (refresh polyfun cols) -> 75 (causal)
#        -> 48 (cross-ancestry) -> 217 (stratified atlas) -> 46d (convergence) -> 27b (presets)
set -eo pipefail
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "host=$(hostname) START=$(date)"

step () { echo ""; echo "================ $1  ($(date +%H:%M:%S)) ================"; }

step "27a assemble_evidence_atlas";   Rscript 27a_assemble_evidence_atlas.R
step "77b add_polyfun_atlas_columns"; Rscript 77b_add_polyfun_atlas_columns.R
step "75 integrate_causal_overhaul";  Rscript 75_integrate_causal_overhaul.R
step "48 cross_ancestry_replication"; Rscript 48_cross_ancestry_replication.R
step "217 stratified_causal_atlas";   Rscript 217_stratified_causal_atlas.R
step "46d convergence_evidence";      Rscript 46d_convergence_evidence.R
step "27b benchmark_presets";         Rscript 27b_benchmark_presets.R

echo ""
echo "CHAIN COMPLETE rc=$? END=$(date)"
