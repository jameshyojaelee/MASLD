#!/bin/bash
#SBATCH --job-name=download_hmsma
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=../logs/download_hmsma_%j.out
#SBATCH --error=../logs/download_hmsma_%j.err

##############################################################################
# 00g: Download HMSMA Processed Visium Data (Li et al., Nat Genet 2025)
#
# STATUS: HMSMA data is ACCESS-CONTROLLED.
#   - Contact: jin.chai@cldcsw.org
#   - Portal: https://db.genomics.cn/stomics/hmsma/download
#   - Embargo note: Some data embargoed until June 30, 2026
#
# This script handles data once access is granted. Two modes:
#   Mode 1: Manual download — user downloads files via browser, places in
#           DATA_DIR, then runs this script with --organize-only
#   Mode 2: Automated download — once download URLs are known, script
#           downloads all files via wget/aria2c
#
# Files needed (175 Stomics + 35 H&E + 1 metadata + 1 scRNA):
#   - 35 × filtered_feature_bc_matrix.h5  (Visium count matrices)
#   - 35 × scalefactors_json.json         (scale factors)
#   - 35 × tissue_lowres_image.png        (low-res tissue image)
#   - 35 × tissue_hires_image.png         (high-res tissue image)
#   - 35 × tissue_positions.csv           (spot coordinates)
#   - 35 × H&E images (.tif or .jpg)      (original staining)
#   - Supplementary_Table_1.docx           (sample metadata)
#   - adata_29samples.h5ad                 (merged scRNA-seq reference)
##############################################################################

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
DATA_DIR="${PROJECT_ROOT}/data/HRA007511/hmsma_processed"
SCRNA_DIR="${PROJECT_ROOT}/data/HRA007511/hmsma_scrna"
METADATA_DIR="${PROJECT_ROOT}/data/HRA007511/metadata"
SPATIAL_ROOT="${PROJECT_ROOT}/Analysis/Spatial"
LOG_DIR="${SPATIAL_ROOT}/logs"

mkdir -p "${DATA_DIR}/stomics"
mkdir -p "${DATA_DIR}/he_images"
mkdir -p "${SCRNA_DIR}"
mkdir -p "${METADATA_DIR}"
mkdir -p "${LOG_DIR}"

# ─── Sample Registry ────────────────────────────────────────────────────────
# 35 Visium samples from HMSMA portal
# Format: CONDITION-PATIENT_ID
CTRL_SAMPLES=(CTRL-161 CTRL-179 CTRL-180 CTRL-5113 CTRL-5759 CTRL-8715 CTRL-8395)
MASLD_SAMPLES=(MASLD-0966 MASLD-1479 MASLD-1492 MASLD-1493 MASLD-1495 MASLD-1497 MASLD-1498 MASLD-2768 MASLD-4973 MASLD-9993)
MASH_SAMPLES=(MASH-0413 MASH-0422 MASH-0835 MASH-1086 MASH-1475 MASH-1478 MASH-1480 MASH-1481 MASH-1494 MASH-1501 MASH-2534 MASH-3096 MASH-3344 MASH-4426 MASH-7866 MASH-8684 MASH-9136 MASH-9440)

ALL_SAMPLES=("${CTRL_SAMPLES[@]}" "${MASLD_SAMPLES[@]}" "${MASH_SAMPLES[@]}")

# Per-sample Stomics files (5 per sample = 175 total)
STOMICS_SUFFIXES=(
    "_filtered_feature_bc_matrix.h5"
    "_scalefactors_json.json"
    "_tissue_lowres_image.png"
    "_tissue_hires_image.png"
    "_tissue_positions.csv"
)

# H&E images — mixed formats (.tif for batch 1, .jpg for batch 2)
# Batch 1 (older, .tif): CTRL-161,179,180 + MASLD/MASL-1479,1492,1493,1495,1497,1498
#   + MASH-1475,1478,1480,1481,1494,1501
# Batch 2 (newer, .jpg): all others
# NOTE: H&E naming uses MASL (not MASLD) for steatosis samples
declare -A HE_EXTENSIONS
for s in CTRL-161 CTRL-179 CTRL-180 MASH-1475 MASH-1478 MASH-1480 MASH-1481 MASH-1494 MASH-1501; do
    HE_EXTENSIONS[$s]=".tif"
done
# MASLD H&E images use MASL prefix (no D)
for s in MASLD-1479 MASLD-1492 MASLD-1493 MASLD-1495 MASLD-1497 MASLD-1498; do
    HE_EXTENSIONS[$s]=".tif"
done
for s in CTRL-5113 CTRL-5759 CTRL-8395 CTRL-8715 \
         MASLD-0966 MASLD-2768 MASLD-4973 MASLD-9993 \
         MASH-0413 MASH-0422 MASH-0835 MASH-1086 MASH-2534 MASH-3096 \
         MASH-3344 MASH-4426 MASH-7866 MASH-8684 MASH-9136 MASH-9440; do
    HE_EXTENSIONS[$s]=".jpg"
done

# ─── Download URL Configuration ─────────────────────────────────────────────
# TODO: Update BASE_URL once HMSMA grants access
# Possible patterns (to be confirmed after access):
#   https://ftp.cngb.org/pub/stomics/STTXXXXXXX/stomics/{filename}
#   Direct links provided by jin.chai@cldcsw.org
BASE_URL="${HMSMA_BASE_URL:-}"  # Set via environment variable

# ─── Functions ──────────────────────────────────────────────────────────────

download_file() {
    local url="$1"
    local output="$2"
    local max_retries=5

    if [[ -f "${output}" ]]; then
        echo "  EXISTS: $(basename ${output})"
        return 0
    fi

    for attempt in $(seq 1 ${max_retries}); do
        echo "  Downloading $(basename ${output}) (attempt ${attempt}/${max_retries})..."
        if wget -q --timeout=120 --tries=1 -c -O "${output}.tmp" "${url}"; then
            mv "${output}.tmp" "${output}"
            echo "  OK: $(basename ${output}) ($(du -h ${output} | cut -f1))"
            return 0
        fi
        sleep $((attempt * 10))
    done

    echo "  FAILED: $(basename ${output}) after ${max_retries} attempts"
    rm -f "${output}.tmp"
    return 1
}

generate_file_list() {
    # Generate list of all files to download (for manual reference)
    local output="${DATA_DIR}/file_list.txt"
    echo "# HMSMA file list — $(date)" > "${output}"
    echo "# 175 Stomics + 35 H&E + metadata + scRNA" >> "${output}"
    echo "" >> "${output}"

    echo "## Stomics (175 files)" >> "${output}"
    for sample in "${ALL_SAMPLES[@]}"; do
        for suffix in "${STOMICS_SUFFIXES[@]}"; do
            echo "${sample}${suffix}" >> "${output}"
        done
    done

    echo "" >> "${output}"
    echo "## H&E Staining (35 files)" >> "${output}"
    for sample in "${ALL_SAMPLES[@]}"; do
        local ext="${HE_EXTENSIONS[$sample]:-.jpg}"
        # MASLD samples use MASL prefix for H&E
        local he_name="${sample}"
        if [[ "${sample}" == MASLD-* ]]; then
            he_name="MASL-${sample#MASLD-}"
        fi
        echo "${he_name}${ext}" >> "${output}"
    done

    echo "" >> "${output}"
    echo "## Metadata" >> "${output}"
    echo "Supplementary_Table_1.docx" >> "${output}"

    echo "" >> "${output}"
    echo "## scRNA-seq (merged)" >> "${output}"
    echo "adata_29samples.h5ad" >> "${output}"

    echo "File list written to: ${output}"
    echo "Total files: $(grep -c '^[^#]' ${output} | head -1)"
}

validate_downloads() {
    local n_ok=0
    local n_missing=0
    local missing_files=()

    echo ""
    echo "=== Validation ==="

    for sample in "${ALL_SAMPLES[@]}"; do
        local all_present=true
        for suffix in "${STOMICS_SUFFIXES[@]}"; do
            local f="${DATA_DIR}/stomics/${sample}${suffix}"
            if [[ ! -f "$f" ]]; then
                all_present=false
                missing_files+=("${sample}${suffix}")
            fi
        done
        if $all_present; then
            local h5_size=$(du -h "${DATA_DIR}/stomics/${sample}_filtered_feature_bc_matrix.h5" 2>/dev/null | cut -f1)
            echo "  ${sample}: OK (${h5_size})"
            ((n_ok++))
        else
            echo "  ${sample}: INCOMPLETE"
            ((n_missing++))
        fi
    done

    echo ""
    echo "Stomics: ${n_ok}/${#ALL_SAMPLES[@]} samples complete"

    # Check H&E
    local n_he=0
    for sample in "${ALL_SAMPLES[@]}"; do
        local ext="${HE_EXTENSIONS[$sample]:-.jpg}"
        local he_name="${sample}"
        if [[ "${sample}" == MASLD-* ]]; then
            he_name="MASL-${sample#MASLD-}"
        fi
        if [[ -f "${DATA_DIR}/he_images/${he_name}${ext}" ]]; then
            ((n_he++))
        fi
    done
    echo "H&E images: ${n_he}/35"

    # Check metadata
    [[ -f "${METADATA_DIR}/Supplementary_Table_1.docx" ]] && echo "Metadata: OK" || echo "Metadata: MISSING"

    # Check scRNA
    [[ -f "${SCRNA_DIR}/adata_29samples.h5ad" ]] && echo "scRNA-seq: OK" || echo "scRNA-seq: MISSING"

    echo ""
    if [[ ${n_missing} -gt 0 ]]; then
        echo "Missing files (${#missing_files[@]}):"
        for f in "${missing_files[@]}"; do
            echo "  ${f}"
        done
    fi

    return ${n_missing}
}

# ─── Main ───────────────────────────────────────────────────────────────────

echo "================================================================="
echo "00g: Download HMSMA Processed Visium Data"
echo "  Samples: ${#ALL_SAMPLES[@]} (${#CTRL_SAMPLES[@]} CTRL, ${#MASLD_SAMPLES[@]} MASLD, ${#MASH_SAMPLES[@]} MASH)"
echo "  Date: $(date)"
echo "================================================================="

# Always generate file list for reference
generate_file_list

# Parse mode
MODE="${1:-validate}"

case "${MODE}" in
    validate|--validate)
        echo ""
        echo "Mode: validate (checking existing files)"
        validate_downloads
        ;;

    download|--download)
        if [[ -z "${BASE_URL}" ]]; then
            echo ""
            echo "ERROR: HMSMA_BASE_URL not set."
            echo ""
            echo "HMSMA data requires access authorization."
            echo "Please contact: jin.chai@cldcsw.org"
            echo "Portal: https://db.genomics.cn/stomics/hmsma/download"
            echo ""
            echo "Once access is granted and you have the download URL, run:"
            echo "  HMSMA_BASE_URL=https://... sbatch 00g_download_hmsma_processed.sh download"
            echo ""
            echo "Alternatively, download files manually and place them in:"
            echo "  Stomics files → ${DATA_DIR}/stomics/"
            echo "  H&E images    → ${DATA_DIR}/he_images/"
            echo "  Metadata      → ${METADATA_DIR}/"
            echo "  scRNA-seq     → ${SCRNA_DIR}/"
            echo "Then run: bash 00g_download_hmsma_processed.sh validate"
            exit 1
        fi

        echo ""
        echo "Mode: download (BASE_URL=${BASE_URL})"

        # Download Stomics files
        echo ""
        echo "--- Stomics (175 files) ---"
        n_fail=0
        for sample in "${ALL_SAMPLES[@]}"; do
            for suffix in "${STOMICS_SUFFIXES[@]}"; do
                fname="${sample}${suffix}"
                download_file "${BASE_URL}/stomics/${fname}" "${DATA_DIR}/stomics/${fname}" || ((n_fail++))
            done
        done

        # Download H&E images
        echo ""
        echo "--- H&E Staining (35 files) ---"
        for sample in "${ALL_SAMPLES[@]}"; do
            ext="${HE_EXTENSIONS[$sample]:-.jpg}"
            he_name="${sample}"
            if [[ "${sample}" == MASLD-* ]]; then
                he_name="MASL-${sample#MASLD-}"
            fi
            fname="${he_name}${ext}"
            download_file "${BASE_URL}/he/${fname}" "${DATA_DIR}/he_images/${fname}" || ((n_fail++))
        done

        # Download metadata
        echo ""
        echo "--- Metadata ---"
        download_file "${BASE_URL}/metadata/Supplementary_Table_1.docx" "${METADATA_DIR}/Supplementary_Table_1.docx" || ((n_fail++))

        # Download scRNA-seq reference
        echo ""
        echo "--- scRNA-seq ---"
        download_file "${BASE_URL}/scrna/adata_29samples.h5ad" "${SCRNA_DIR}/adata_29samples.h5ad" || ((n_fail++))

        echo ""
        echo "Download complete. Failures: ${n_fail}"
        validate_downloads
        ;;

    file-list|--file-list)
        echo "File list generated at: ${DATA_DIR}/file_list.txt"
        ;;

    *)
        echo "Usage: $0 {validate|download|file-list}"
        echo ""
        echo "  validate   — Check which files are present (default)"
        echo "  download   — Download files (requires HMSMA_BASE_URL env var)"
        echo "  file-list  — Generate file list only"
        exit 1
        ;;
esac

echo ""
echo "Done: $(date)"
