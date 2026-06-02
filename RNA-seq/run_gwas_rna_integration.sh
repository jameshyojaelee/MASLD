#!/bin/bash
# run_gwas_rna_integration.sh — SLURM orchestrator for GWAS-RNA-scRNA co-analysis
#
# Phase A: Computational analyses (200-205)
# Phase B: Integration (206) — depends on all Phase A
#
# Usage: bash run_gwas_rna_integration.sh
#
# All scripts output to: RNA-seq/results/gwas_rna_integration/

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
LOGDIR="${BASE}/RNA-seq/logs/gwas_rna_integration"
mkdir -p "$LOGDIR"

export MASLD_PROJECT_ROOT="$BASE"

echo "=== GWAS-RNA-scRNA Co-Analysis Pipeline ==="
echo "Start: $(date)"
echo "Logs: $LOGDIR"
echo ""

# ---------- Script 200: INTACT Composite Scoring ----------
JOB_200=$(sbatch --parsable \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=32G \
  --time=48:00:00 \
  --job-name=intact_200 \
  --output="${LOGDIR}/200_intact_%j.out" \
  --error="${LOGDIR}/200_intact_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${BASE}/RNA-seq/200_intact_scoring.R")
echo "Submitted Script 200 (INTACT): Job $JOB_200"

# ---------- Script 201: Cell-Type Heritability Enrichment ----------
JOB_201=$(sbatch --parsable \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=64G \
  --time=48:00:00 \
  --job-name=celltype_herit_201 \
  --output="${LOGDIR}/201_celltype_herit_%j.out" \
  --error="${LOGDIR}/201_celltype_herit_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${BASE}/RNA-seq/201_celltype_heritability.R")
echo "Submitted Script 201 (Cell-type heritability): Job $JOB_201"

# ---------- Script 202: Gene-Set Enrichment (depends on 201 for gene sets) ----------
JOB_202=$(sbatch --parsable \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=32G \
  --time=48:00:00 \
  --job-name=geneset_enrich_202 \
  --output="${LOGDIR}/202_geneset_enrich_%j.out" \
  --error="${LOGDIR}/202_geneset_enrich_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${BASE}/RNA-seq/202_geneset_enrichment.R")
echo "Submitted Script 202 (Gene-set enrichment): Job $JOB_202"

# ---------- Script 203: (Removed — MR excluded from pipeline) ----------

# ---------- Script 204: scDRS Scoring ----------
JOB_204=$(sbatch --parsable \
  --partition=gpu \
  --gres=gpu:1 \
  --cpus-per-task=8 \
  --mem=64G \
  --time=48:00:00 \
  --job-name=scdrs_204 \
  --output="${LOGDIR}/204_scdrs_%j.out" \
  --error="${LOGDIR}/204_scdrs_%j.err" \
  --wrap="micromamba run -n rapids_singlecell python ${BASE}/RNA-seq/204_scdrs_scoring.py")
echo "Submitted Script 204 (scDRS): Job $JOB_204"

# ---------- Script 205: Proportion × COLOC ----------
JOB_205=$(sbatch --parsable \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=32G \
  --time=48:00:00 \
  --job-name=prop_coloc_205 \
  --output="${LOGDIR}/205_prop_coloc_%j.out" \
  --error="${LOGDIR}/205_prop_coloc_%j.err" \
  --wrap="micromamba run -n rnaseq Rscript ${BASE}/RNA-seq/205_proportion_coloc.R")
echo "Submitted Script 205 (Proportion × COLOC): Job $JOB_205"

echo ""
echo "=== All Phase A jobs submitted ==="
echo "Jobs: $JOB_200, $JOB_201, $JOB_202, $JOB_203, $JOB_204, $JOB_205"
echo ""
echo "Monitor with: squeue -u \$USER --name=intact_200,celltype_herit_201,geneset_enrich_202,mr_mediation_203,scdrs_204,prop_coloc_205"
echo ""
echo "After all complete, run Phase B integration script manually:"
echo "  Rscript ${BASE}/RNA-seq/206_gwas_rna_integration.R"
