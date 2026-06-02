#!/bin/bash
#SBATCH --job-name=dl_liver_gwas
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/download_liver_gwas_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/download_liver_gwas_%j.err

# Download UKBB liver enzyme GWAS summary statistics from EBI GWAS Catalog
# AST (GCST90019497), GGT (GCST90019507) — same study as ALT (Sinnott-Armstrong 2021)
# PDFF (GCST90267352) — Pazoki 2022, GRCh37

set -euo pipefail

OUTDIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data"
mkdir -p "$OUTDIR"

echo "=== Downloading UKBB liver enzyme GWAS summary statistics ==="
echo "Start: $(date)"

# ---------------------------------------------------------------------------
# AST (GCST90019497) — Sinnott-Armstrong 2021, N=343,850, hg38 harmonised
# ---------------------------------------------------------------------------
AST_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90019001-GCST90020000/GCST90019497/harmonised/GCST90019497.h.tsv.gz"
AST_OUT="$OUTDIR/GCST90019497_UKBB_AST_harmonised.tsv.gz"

if [ -f "$AST_OUT" ]; then
  echo "AST already downloaded: $AST_OUT"
else
  echo "Downloading AST (GCST90019497)..."
  wget -c -q --show-progress -O "$AST_OUT" "$AST_URL"
  echo "  Saved: $AST_OUT ($(du -h "$AST_OUT" | cut -f1))"
fi

# ---------------------------------------------------------------------------
# GGT (GCST90019507) — Sinnott-Armstrong 2021, N=343,850, hg38 harmonised
# ---------------------------------------------------------------------------
GGT_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90019001-GCST90020000/GCST90019507/harmonised/GCST90019507.h.tsv.gz"
GGT_OUT="$OUTDIR/GCST90019507_UKBB_GGT_harmonised.tsv.gz"

if [ -f "$GGT_OUT" ]; then
  echo "GGT already downloaded: $GGT_OUT"
else
  echo "Downloading GGT (GCST90019507)..."
  wget -c -q --show-progress -O "$GGT_OUT" "$GGT_URL"
  echo "  Saved: $GGT_OUT ($(du -h "$GGT_OUT" | cut -f1))"
fi

# ---------------------------------------------------------------------------
# PDFF (GCST90267352) — Pazoki 2022, N=33,588, GRCh37
# Note: This uses GWAS-SSF format and is GRCh37 (hg19), not hg38
# ---------------------------------------------------------------------------
PDFF_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267352/GCST90267352.tsv.gz"
PDFF_OUT="$OUTDIR/GCST90267352_PDFF_Pazoki2022.tsv.gz"

if [ -f "$PDFF_OUT" ]; then
  echo "PDFF already downloaded: $PDFF_OUT"
else
  echo "Downloading PDFF (GCST90267352)..."
  wget -c -q --show-progress -O "$PDFF_OUT" "$PDFF_URL"
  echo "  Saved: $PDFF_OUT ($(du -h "$PDFF_OUT" | cut -f1))"
fi

# ---------------------------------------------------------------------------
# Verify downloads
# ---------------------------------------------------------------------------
echo ""
echo "=== Verifying downloads ==="
for f in "$AST_OUT" "$GGT_OUT" "$PDFF_OUT"; do
  if [ -f "$f" ]; then
    SIZE=$(du -h "$f" | cut -f1)
    # Check gzip integrity
    if gzip -t "$f" 2>/dev/null; then
      echo "  OK: $(basename "$f") ($SIZE)"
    else
      echo "  WARN: $(basename "$f") ($SIZE) — gzip integrity check failed"
    fi
  else
    echo "  MISSING: $(basename "$f")"
  fi
done

# Quick peek at column headers
echo ""
echo "=== Column headers ==="
echo "AST:"
zcat "$AST_OUT" 2>/dev/null | head -1 || echo "  (could not read)"
echo "GGT:"
zcat "$GGT_OUT" 2>/dev/null | head -1 || echo "  (could not read)"
echo "PDFF:"
zcat "$PDFF_OUT" 2>/dev/null | head -1 || echo "  (could not read)"

echo ""
echo "=== Download complete ==="
echo "End: $(date)"
