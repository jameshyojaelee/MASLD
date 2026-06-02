#!/bin/bash
#SBATCH --job-name=bbj_download
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/BBJ/logs/download_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/BBJ/logs/download_%j.err

# Download Biobank Japan (BBJ) Liver Enzyme GWAS Summary Statistics
#
# Source: PheWeb Japan (https://pheweb.jp/)
# Publication: Kanai et al. 2018, Nature Genetics (PMID: 30108127)
# Sample sizes: N ~ 160,000-260,000 (Japanese ancestry)
#
# These are quantitative trait GWAS for liver enzymes:
#   - ALT (alanine aminotransferase)
#   - AST (aspartate aminotransferase)
#   - GGT (gamma-glutamyl transferase)
#
# Genome build: GRCh37 (hg19) — requires liftover to GRCh38
# Format: PheWeb standard (tab-separated)
#
# After download, run format_bbj_for_coloc.R to harmonize columns and liftover.

set -euo pipefail

BBJ_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/BBJ"
mkdir -p "$BBJ_DIR/logs" "$BBJ_DIR/raw"

echo "============================================================"
echo "  BBJ Liver Enzyme GWAS Download"
echo "  Started: $(date)"
echo "  Output: $BBJ_DIR/raw/"
echo "============================================================"
echo ""

# PheWeb Japan download URLs (Kanai et al. 2018)
# These are publicly available without registration
# URL pattern: https://pheweb.jp/download/[pheno]
#
# NOTE: The exact download mechanism for PheWeb Japan may require
# checking the current website. Below are the expected URLs based
# on the PheWeb Japan API. If these fail, visit https://pheweb.jp/
# and navigate to the phenotype pages for ALT, AST, GGT to find
# the current download links.

declare -A PHENOS
PHENOS[ALT]="ALT"
PHENOS[AST]="AST"
PHENOS[GGT]="GGT"

# BBJ sample sizes from Kanai et al. 2018, Nature Genetics
declare -A SAMPLE_SIZES
SAMPLE_SIZES[ALT]=261406
SAMPLE_SIZES[AST]=261270
SAMPLE_SIZES[GGT]=164006

for TRAIT in ALT AST GGT; do
    PHENO="${PHENOS[$TRAIT]}"
    OUTFILE="$BBJ_DIR/raw/BBJ_${TRAIT}_raw.tsv.gz"

    if [ -f "$OUTFILE" ]; then
        echo "  $TRAIT: already downloaded ($OUTFILE)"
        continue
    fi

    echo "  Downloading BBJ $TRAIT (N=${SAMPLE_SIZES[$TRAIT]})..."

    # Primary: PheWeb Japan download endpoint
    # The BBJ GWAS are also available via GWAS Catalog or JENGER
    # Try PheWeb Japan first, then fall back to alternative sources
    URL="https://pheweb.jp/download/${PHENO}"

    if wget -q --spider "$URL" 2>/dev/null; then
        wget -O "$OUTFILE" "$URL"
        echo "    Downloaded from PheWeb Japan: $(ls -lh "$OUTFILE" | awk '{print $5}')"
    else
        echo "    WARNING: PheWeb Japan URL may have changed for $TRAIT"
        echo "    Try these alternative sources:"
        echo "      1. https://pheweb.jp/ — navigate to $TRAIT phenotype"
        echo "      2. JENGER: http://jenger.riken.jp/en/result — BBJ quantitative traits"
        echo "      3. GWAS Catalog: search 'Kanai 2018' + $TRAIT"
        echo ""
        echo "    After manual download, place file at: $OUTFILE"
        echo "    Expected format: tab-separated with columns including"
        echo "    chrom, pos, ref, alt, pval, beta, sebeta (or similar)"
    fi
done

echo ""
echo "  Saving sample size metadata..."
cat > "$BBJ_DIR/bbj_sample_sizes.tsv" << 'EOF'
trait	N	ancestry	publication	pmid
ALT	261406	East_Asian	Kanai et al. 2018 Nat Genet	30108127
AST	261270	East_Asian	Kanai et al. 2018 Nat Genet	30108127
GGT	164006	East_Asian	Kanai et al. 2018 Nat Genet	30108127
EOF

echo ""
echo "============================================================"
echo "  Download complete: $(date)"
echo "  Next step: Rscript format_bbj_for_coloc.R"
echo "============================================================"
