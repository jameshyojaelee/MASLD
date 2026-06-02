#!/bin/bash
#SBATCH --partition=bigmem --qos=interactive
#SBATCH --cpus-per-task=16 --mem=200G --time=48:00:00
#SBATCH --job-name=nmf95_nonprotonly
#SBATCH --output=RNA-seq/logs/nmf95_nonprotonly_%j.out
#SBATCH --error=RNA-seq/logs/nmf95_nonprotonly_%j.err
# Bravo Task 7 — PROT_ONLY=FALSE control refit.
# Extended strip (MT-pseudogenes + expanded X-escape + HLA + pseudogene biotype
# + UTX/KDM6A dup fix) WITHOUT the protein-coding-only gate.
# K_SWEEP {3..8} fit inline (one bigmem job, ~3-5 hr expected at nrun=20).
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export NMF_PROTEIN_CODING_ONLY=FALSE
export NMF_STRIP_HLA=TRUE
export NMF_INLINE_ALL_K=TRUE
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean_nonprotonly.rds"

echo "=== PROT_ONLY=FALSE remediation refit ==="
echo "  Cache  : $NMF_CACHE_PATH"
echo "  Strip  : HLA=TRUE, PROT_ONLY=FALSE, pseudogenes+MT-pseudo+P\\d+-suffix"
echo "  K_SWEEP: 3..8 inline (no shards)"
Rscript RNA-seq/95_nmf_clean_ksweep.R
