#!/bin/bash
#SBATCH --job-name=atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/atlas_rebuild_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/atlas_rebuild_%j.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=48:00:00
# Atlas + convergence rebuild after the 46d genetic-gate fix (statistical audit 2026-06-29):
#   - 46d S2a gate now built ONLY from canonical SuSiE-with-ABF-fallback (coloc_susie_best_pp4 /
#     coloc_abf_best_pp4); legacy sc-eQTL / per-trait _coloc_pp4 columns removed.
#   - sex_interaction modality dropped from the S1 disease-direction consensus.
# Chain: 27a (assemble atlas) -> 75 (causal integrate) -> 217 (stratified causal) -> 46d (convergence).
set -eo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "$MASLD_PROJECT_ROOT/RNA-seq/logs"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-8}
cd "$MASLD_PROJECT_ROOT/RNA-seq"
echo "[$(date)] START 27a"; Rscript 27a_assemble_evidence_atlas.R
echo "[$(date)] START 75";  Rscript 75_integrate_causal_overhaul.R
echo "[$(date)] START 217"; Rscript 217_stratified_causal_atlas.R
echo "[$(date)] START 46d"; Rscript 46d_convergence_evidence.R
echo "[$(date)] DONE atlas+convergence rebuild"
