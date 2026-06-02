#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=gwas_atac_56_motif
#SBATCH --output=../logs/56_motif_%j.out
#SBATCH --error=../logs/56_motif_%j.err

source ~/.bashrc
micromamba activate motifbreakr
module load meme/5.5.9

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src"
cd "$SCRIPT_DIR"

# FIX (review Script-56 drift): run v2 (canonical — re-resolves REF/ALT from BSgenome, the F100
# dependency for 56c). MB_WRITE_CANONICAL=1 makes v2 also write the unsuffixed
# motif_disruption_scores.csv that 56c/57/56i consume. Background defaults to uniform; set
# MB_BACKGROUND=peak (after running 56e_gc_background.R) for the GC/peak-matched canonical run.
echo "Running 56_motif_disruption_v2.R (canonical BSgenome REF/ALT producer) on $(hostname) at $(date)"
MB_WRITE_CANONICAL=1 Rscript 56_motif_disruption_v2.R
echo "Done at $(date)"
