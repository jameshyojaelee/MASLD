#!/usr/bin/env bash
#SBATCH --job-name=gsmap_slice_mean
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/gsmap_slice_mean_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/gsmap_slice_mean_%j.err
##############################################################################
# 15d_create_slice_mean.sh — Create gsMap slice means for each dataset.
#
# gsmap create_slice_mean requires:
#   --sample_name_list (space-separated sample names)
#   --h5ad_yaml OR --h5ad_list (sample→path mapping)
#   --slice_mean_output_file
#   --data_layer
##############################################################################
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate gsmap

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
INPUT_DIR="${PROJECT_ROOT}/Analysis/Spatial/data/gsmap_input"

cd "${PROJECT_ROOT}"

echo "======================================================================="
echo "  15d: Create gsMap Slice Means"
echo "======================================================================="

# ── GSE192741 ─────────────────────────────────────────────────────────────
GSE_YAML="${INPUT_DIR}/gse192741_h5ad.yaml"
GSE_OUT="${INPUT_DIR}/gse192741_slice_mean.parquet"

if [ -f "${GSE_YAML}" ]; then
    if [ -f "${GSE_OUT}" ]; then
        echo "  GSE192741 slice mean already exists: ${GSE_OUT}"
    else
        echo ""
        echo "  Creating slice mean for GSE192741..."
        # Extract sample names from YAML (keys before the colon)
        GSE_SAMPLES=$(grep -oP '^\S+(?=:)' "${GSE_YAML}" | tr '\n' ' ')
        echo "  Samples: ${GSE_SAMPLES}"
        gsmap create_slice_mean \
            --sample_name_list ${GSE_SAMPLES} \
            --h5ad_yaml "${GSE_YAML}" \
            --slice_mean_output_file "${GSE_OUT}" \
            --data_layer counts
        echo "  GSE192741 slice mean saved: ${GSE_OUT}"
    fi
else
    echo "  WARNING: ${GSE_YAML} not found — run 15c first"
fi

# ── Vu et al. ─────────────────────────────────────────────────────────────
VU_YAML="${INPUT_DIR}/vu_h5ad.yaml"
VU_OUT="${INPUT_DIR}/vu_slice_mean.parquet"

if [ -f "${VU_YAML}" ]; then
    if [ -f "${VU_OUT}" ]; then
        echo "  Vu slice mean already exists: ${VU_OUT}"
    else
        echo ""
        echo "  Creating slice mean for Vu et al...."
        VU_SAMPLES=$(grep -oP '^\S+(?=:)' "${VU_YAML}" | tr '\n' ' ')
        echo "  Samples: ${VU_SAMPLES}"
        gsmap create_slice_mean \
            --sample_name_list ${VU_SAMPLES} \
            --h5ad_yaml "${VU_YAML}" \
            --slice_mean_output_file "${VU_OUT}" \
            --data_layer counts
        echo "  Vu slice mean saved: ${VU_OUT}"
    fi
else
    echo "  WARNING: ${VU_YAML} not found — run 15c first"
fi

echo ""
echo "======================================================================="
echo "  15d: Complete"
echo "======================================================================="
echo ""
echo "  Outputs:"
ls -lh "${INPUT_DIR}"/*slice_mean* 2>/dev/null || echo "  (none)"
