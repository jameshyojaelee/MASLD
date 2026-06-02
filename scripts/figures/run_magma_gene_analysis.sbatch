#!/bin/bash
#SBATCH --job-name=magma_ga
#SBATCH --partition=cpu
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --array=0-18
#SBATCH --output=scripts/figures/logs/magma_ga_%A_%a.out
#SBATCH --error=scripts/figures/logs/magma_ga_%A_%a.err

# Run MAGMA gene-analysis (one GWAS per SLURM array task).
# Input list comes from tools/magma/work/gwas_inputs/_N_table.tsv built by
# run_magma_pipeline.sbatch. Each task reads its (gwas_id, N) row by index.

set -eo pipefail
mkdir -p scripts/figures/logs

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

MAGMA_BIN=tools/magma/magma
WORK=tools/magma/work
REF=tools/magma/ref
INPUT_DIR=$WORK/gwas_inputs
OUT_DIR=$WORK/gene_analysis_out
mkdir -p $OUT_DIR

# Read the (gwas_id, N, status) for this task from the audit table.
# IMPORTANT: use tab as separator — the "family" field (e.g. "Liver enzymes
# (EUR)") contains spaces and breaks default-whitespace awk splitting.
N_TABLE=$INPUT_DIR/_N_table.tsv
ROW=$(awk -F'\t' -v idx=$SLURM_ARRAY_TASK_ID 'NR==idx+2 {print $1"\t"$3"\t"$8}' $N_TABLE)
GWAS_ID=$(echo "$ROW" | cut -f1)
N_VALUE=$(echo "$ROW" | cut -f2)
STATUS=$(echo "$ROW"  | cut -f3)

echo "[gene_analysis $SLURM_ARRAY_TASK_ID] gwas_id=$GWAS_ID  N=$N_VALUE  status=$STATUS"

if [ "$STATUS" != "OK" ]; then
  echo "  SKIP — non-OK status"
  exit 0
fi

INPUT=$INPUT_DIR/${GWAS_ID}.magma_input.tsv
if [ ! -s "$INPUT" ]; then
  echo "  SKIP — missing or empty input: $INPUT"
  exit 0
fi

# Run MAGMA gene-analysis using the pre-built annotation + 1000G LD bfile
$MAGMA_BIN \
  --bfile $REF/g1000_eur \
  --pval $INPUT use=SNP,P ncol=N \
  --gene-annot $WORK/magma_annot.genes.annot \
  --out $OUT_DIR/$GWAS_ID

echo "[gene_analysis $SLURM_ARRAY_TASK_ID] done; output: $OUT_DIR/${GWAS_ID}.genes.out"
