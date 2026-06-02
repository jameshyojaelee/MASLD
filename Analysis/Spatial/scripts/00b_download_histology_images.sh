#!/usr/bin/env bash
#SBATCH --job-name=dl_histology
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/download_histology_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/download_histology_%j.err
##############################################################################
# 00b_download_histology_images.sh — Download H&E tissue images
#
# GSE192741: Extract tissue images from GEO GSE192741_RAW.tar
#   - Images bundled as tissue_hires_image_{SAMPLE}.png.gz per sample
#   - Already present in SpaceRanger outputs; this extracts standalone copies
#
# HRA007511: NGDC does NOT deposit H&E images (FASTQs only)
#   - SpaceRanger will use --unknown-slide (expression-only tissue detection)
#   - This is standard and does not affect downstream analysis
##############################################################################
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

SPATIAL_ROOT="Analysis/Spatial"
GSE_DIR="${SPATIAL_ROOT}/metadata/histology/GSE192741"
HRA_DIR="${SPATIAL_ROOT}/metadata/histology/HRA007511"
mkdir -p "${GSE_DIR}" "${HRA_DIR}"

echo "=============================================="
echo "  Histology Image Download"
echo "  Date: $(date)"
echo "=============================================="

# ── GSE192741 ──────────────────────────────────────────────────────────────
echo ""
echo "=== GSE192741 (GEO RAW.tar) ==="

# Images are inside GSE192741_RAW.tar as gzipped PNGs:
#   GSM5764414_tissue_hires_image_JBO001.png.gz
#   GSM5764414_tissue_lowres_image_JBO001.png.gz
# etc. for all 15 samples.

# First check if images already exist from SpaceRanger outputs
N_EXISTING=$(find "${SPATIAL_ROOT}/results/spaceranger/GSE192741" \
    -name "tissue_hires_image.png" 2>/dev/null | wc -l)
echo "  Images in SpaceRanger outputs: ${N_EXISTING}"

if [ "${N_EXISTING}" -gt 0 ]; then
    echo "  Images already available in SpaceRanger outs/spatial/ directories."
    echo "  Creating symlinks to standalone histology directory..."
    for sr_dir in ${SPATIAL_ROOT}/results/spaceranger/GSE192741/*/outs/spatial; do
        if [ -f "${sr_dir}/tissue_hires_image.png" ]; then
            SAMPLE=$(basename "$(dirname "$(dirname "${sr_dir}")")")
            ln -sf "$(realpath "${sr_dir}/tissue_hires_image.png")" \
                "${GSE_DIR}/${SAMPLE}_hires.png" 2>/dev/null || true
            ln -sf "$(realpath "${sr_dir}/tissue_lowres_image.png")" \
                "${GSE_DIR}/${SAMPLE}_lowres.png" 2>/dev/null || true
        fi
    done
else
    echo "  No SpaceRanger outputs found. Downloading from GEO..."
    GEO_TAR_URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE192nnn/GSE192741/suppl/GSE192741_RAW.tar"
    TMPDIR_TAR=$(mktemp -d)

    echo "  Downloading GSE192741_RAW.tar (~180 MB)..."
    wget -q -O "${TMPDIR_TAR}/GSE192741_RAW.tar" "${GEO_TAR_URL}" 2>&1 || {
        echo "  ERROR: Failed to download GSE192741_RAW.tar"
        rm -rf "${TMPDIR_TAR}"
    }

    if [ -f "${TMPDIR_TAR}/GSE192741_RAW.tar" ]; then
        echo "  Extracting tissue images..."
        # Extract only image files (not h5 matrices)
        tar -xf "${TMPDIR_TAR}/GSE192741_RAW.tar" -C "${TMPDIR_TAR}" \
            --wildcards '*tissue_hires_image*' '*tissue_lowres_image*' 2>/dev/null || true

        # Decompress and rename to sample-level files
        for img_gz in "${TMPDIR_TAR}"/*tissue_hires_image*.png.gz; do
            [ -f "${img_gz}" ] || continue
            # Extract sample name: GSM5764414_tissue_hires_image_JBO001.png.gz → JBO001
            BASENAME=$(basename "${img_gz}" .png.gz)
            SAMPLE=$(echo "${BASENAME}" | sed 's/.*tissue_hires_image_//')
            gunzip -c "${img_gz}" > "${GSE_DIR}/${SAMPLE}_hires.png"
            echo "    Extracted: ${SAMPLE}_hires.png"
        done
        for img_gz in "${TMPDIR_TAR}"/*tissue_lowres_image*.png.gz; do
            [ -f "${img_gz}" ] || continue
            BASENAME=$(basename "${img_gz}" .png.gz)
            SAMPLE=$(echo "${BASENAME}" | sed 's/.*tissue_lowres_image_//')
            gunzip -c "${img_gz}" > "${GSE_DIR}/${SAMPLE}_lowres.png"
        done

        rm -rf "${TMPDIR_TAR}"
    fi
fi

N_GSE=$(find "${GSE_DIR}" -name "*_hires.png" -o -name "*_lowres.png" 2>/dev/null | wc -l)
echo "  GSE192741 image files: ${N_GSE}"

# ── HRA007511 ──────────────────────────────────────────────────────────────
echo ""
echo "=== HRA007511 (NGDC) ==="
echo "  NGDC deposit for HRA007511 contains FASTQs only (verified via FTP listing)."
echo "  No H&E histology images were deposited."
echo ""
echo "  Checking alternative sources..."

# Try NGDC supplementary (unlikely but worth checking)
HRA_FOUND=0
for url_pattern in \
    "https://download.cncb.ac.cn/gsa-human/HRA007511/suppl/" \
    "https://ngdc.cncb.ac.cn/omix/release/HRA007511/" \
    "ftp://download.big.ac.cn/gsa-human/HRA007511/suppl/"; do
    echo "  Trying: ${url_pattern}"
    if wget -q --spider "${url_pattern}" 2>/dev/null; then
        echo "    Found! Downloading image files..."
        wget -q -r -np -nH --cut-dirs=3 \
            -A "*.tif,*.tiff,*.jpg,*.png,*.jpeg" \
            "${url_pattern}" -P "${HRA_DIR}" 2>&1 || true
        HRA_FOUND=1
        break
    else
        echo "    Not available."
    fi
done

N_HRA=$(find "${HRA_DIR}" -name "*.tif" -o -name "*.tiff" -o -name "*.jpg" -o -name "*.png" 2>/dev/null | wc -l)
echo ""
echo "  HRA007511 images found: ${N_HRA}"

if [ "${N_HRA}" -eq 0 ]; then
    echo ""
    echo "  CONCLUSION: No H&E images available for HRA007511."
    echo "  This is expected — NGDC (GSA-Human) deposits typically include only sequencing data."
    echo "  SpaceRanger will use --unknown-slide visium-1 (expression-based tissue detection)."
    echo "  Impact: No H&E overlay in figures, but all computational analysis is unaffected."
    echo ""
    echo "  To obtain images, contact the authors (Li et al., Nat Genet 2025) directly."
fi

# ── Summary ────────────────────────────────────────────────────────────────
echo ""
echo "=============================================="
echo "  Summary"
echo "  GSE192741: ${N_GSE} image files (${N_GSE} / 2 = $((N_GSE / 2)) samples)"
echo "  HRA007511: ${N_HRA} image files"
echo ""
if [ "${N_HRA}" -eq 0 ]; then
    echo "  Pipeline impact:"
    echo "    - GSE192741 (healthy): Full H&E overlay available"
    echo "    - HRA007511 (MASLD):   Expression-only tissue detection (--unknown-slide)"
    echo "    - All downstream analysis (cell2location, zonation, SVGs) unaffected"
fi
echo "  Date: $(date)"
echo "=============================================="
