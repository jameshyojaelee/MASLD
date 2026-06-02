#!/bin/bash
#SBATCH --job-name=magma_pipe
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=scripts/figures/logs/magma_pipe_%j.out
#SBATCH --error=scripts/figures/logs/magma_pipe_%j.err

# One-shot driver for the MAGMA-anchored scDRS prep pipeline:
#   (1) Build SNP loc + chr:pos→rsID lookup from 1000G .bim
#   (2) MAGMA --annotate (SNP→gene mapping; one-time)
#   (3) Prep per-GWAS MAGMA input files with documented N
# This script does NOT run scDRS scoring (separate SLURM array on subsample atlas).

set -eo pipefail
mkdir -p scripts/figures/logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

MAGMA_BIN=tools/magma/magma
WORK=tools/magma/work
REF=tools/magma/ref
mkdir -p $WORK

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

# ─── Step 1: build SNP loc + chr:pos→rsID lookup ───────────────────────────
echo "[pipe] step 1: SNP loc + chr:pos lookup"
python scripts/figures/magma_setup.py

# ─── Step 2: MAGMA --annotate (SNP → gene mapping) ─────────────────────────
echo "[pipe] step 2: MAGMA annotate"
$MAGMA_BIN \
  --annotate \
  --snp-loc $WORK/g1000_eur.snploc \
  --gene-loc $REF/NCBI37.3.gene.loc \
  --out $WORK/magma_annot

# ─── Step 3: prep per-GWAS MAGMA input files ───────────────────────────────
echo "[pipe] step 3: per-GWAS input prep"
python scripts/figures/magma_gwas_prep.py

echo "[pipe] all setup steps complete."
echo
echo "Audit table:"
cat $WORK/gwas_inputs/_N_table.tsv | column -t -s $'\t'
