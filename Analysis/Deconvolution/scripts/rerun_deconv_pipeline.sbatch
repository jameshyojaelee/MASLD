#!/usr/bin/env bash
#SBATCH --job-name=deconv_rerun
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Deconvolution/logs/rerun_deconv_%j.out
#SBATCH --error=Analysis/Deconvolution/logs/rerun_deconv_%j.err

# Re-run full deconvolution pipeline with new ScaleSC reference (1.17M cells)
# Steps:
#   1. Subsample reference to ~100K cells
#   2. Convert to SCE
#   3. Prepare bulk counts for 4 new datasets
#   4. Run MuSiC for all 10 human datasets
#   5. Run BayesPrism for all 10 human datasets

set -euo pipefail

echo "=== Deconvolution Pipeline Re-run ==="
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Date: $(date)"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export R_LIBS_USER=""

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

BASE="Analysis/Deconvolution"
REF_DIR="${BASE}/reference"
SCRIPTS="${BASE}/scripts"
RESULTS="${BASE}/results"
BULK="${BASE}/bulk"

ANNOTATED_H5AD="Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
SUBSAMPLED_H5AD="${REF_DIR}/reference_human_scalesc_100k.h5ad"
NEW_SCE="${REF_DIR}/reference_human_scalesc_sce.rds"

mkdir -p "${REF_DIR}" "${RESULTS}" "${BULK}" "Analysis/Deconvolution/logs"

# =========================================================================
# Step 1: Subsample reference to ~100K cells
# =========================================================================
echo ""
echo "=== Step 1: Subsample reference ==="
if [ -f "${SUBSAMPLED_H5AD}" ]; then
    echo "  Subsampled reference exists, skipping"
else
    # Need rapids_singlecell env for anndata with the annotated h5ad
    micromamba activate rapids_singlecell
    python "${SCRIPTS}/subsample_reference.py" \
        --input "${ANNOTATED_H5AD}" \
        --output "${SUBSAMPLED_H5AD}" \
        --max-cells 100000
    micromamba activate rnaseq
    export R_LIBS_USER=""
fi

# =========================================================================
# Step 2: Convert to SCE
# =========================================================================
echo ""
echo "=== Step 2: Convert to SCE ==="
if [ -f "${NEW_SCE}" ]; then
    echo "  SCE exists, skipping"
else
    Rscript "${SCRIPTS}/05_make_sce.R" \
        "${SUBSAMPLED_H5AD}" \
        "${NEW_SCE}" \
        "cell_type"
fi

# =========================================================================
# Step 3: Prepare bulk counts for 4 new datasets
# =========================================================================
echo ""
echo "=== Step 3: Prepare bulk counts ==="

HUMAN_BASE="RNA-seq/Human/Patient_Cohorts"

declare -A NEW_DATASETS
NEW_DATASETS=(
    ["GSE162694"]="${HUMAN_BASE}/results/GSE162694/counts/featurecounts/gene_counts.txt|${HUMAN_BASE}/pipelines/custom/GSE162694/metadata/samples.tsv"
    ["GSE174478"]="${HUMAN_BASE}/results/GSE174478/counts/featurecounts/gene_counts.txt|${HUMAN_BASE}/pipelines/custom/GSE174478/metadata/samples.tsv"
    ["GSE193066"]="${HUMAN_BASE}/results/GSE193066/counts/featurecounts/gene_counts.txt|${HUMAN_BASE}/pipelines/custom/GSE193066/metadata/samples_with_paths.tsv"
    ["GSE240729"]="${HUMAN_BASE}/results/GSE240729/counts/featurecounts/gene_counts.txt|${HUMAN_BASE}/pipelines/custom/GSE240729/metadata/samples.tsv"
)

for ds in "${!NEW_DATASETS[@]}"; do
    IFS='|' read -r counts_path meta_path <<< "${NEW_DATASETS[$ds]}"
    bulk_dir="${BULK}/${ds}"
    if [ -f "${bulk_dir}/${ds}_counts.tsv" ]; then
        echo "  ${ds}: bulk counts exist, skipping"
    else
        echo "  ${ds}: preparing bulk counts..."
        Rscript "${SCRIPTS}/06_prepare_bulk_counts.R" \
            "${counts_path}" \
            "${meta_path}" \
            "human" \
            "${bulk_dir}" \
            "${ds}"
    fi
done

# =========================================================================
# Step 4: Run MuSiC for all 10 human datasets
# =========================================================================
echo ""
echo "=== Step 4: MuSiC deconvolution ==="

# All 10 human datasets with their bulk count locations
declare -A ALL_DATASETS
ALL_DATASETS=(
    ["GSE126848"]="${BULK}/GSE126848"
    ["GSE130970"]="${BULK}/GSE130970"
    ["GSE135251"]="${BULK}/GSE135251"
    ["GSE167523"]="${BULK}/GSE167523"
    ["PRJNA512027"]="${BULK}/PRJNA512027"
    ["GSE213621"]="${BULK}/GSE213621"
    ["GSE162694"]="${BULK}/GSE162694"
    ["GSE174478"]="${BULK}/GSE174478"
    ["GSE193066"]="${BULK}/GSE193066"
    ["GSE240729"]="${BULK}/GSE240729"
)

for ds in "${!ALL_DATASETS[@]}"; do
    bulk_dir="${ALL_DATASETS[$ds]}"
    counts_tsv="${bulk_dir}/${ds}_counts.tsv"
    meta_tsv="${bulk_dir}/${ds}_metadata.tsv"
    result_dir="${RESULTS}/${ds}"
    mkdir -p "${result_dir}"

    if [ ! -f "${counts_tsv}" ]; then
        echo "  ${ds}: SKIPPING (no bulk counts at ${counts_tsv})"
        continue
    fi

    echo "  ${ds}: running MuSiC..."
    Rscript "${SCRIPTS}/07_run_music.R" \
        "${NEW_SCE}" \
        "${counts_tsv}" \
        "${meta_tsv}" \
        "${result_dir}" \
        "${ds}" || echo "  WARNING: MuSiC failed for ${ds}"
done

# =========================================================================
# Step 5: Run BayesPrism for all 10 human datasets
# =========================================================================
echo ""
echo "=== Step 5: BayesPrism deconvolution ==="

for ds in "${!ALL_DATASETS[@]}"; do
    bulk_dir="${ALL_DATASETS[$ds]}"
    counts_tsv="${bulk_dir}/${ds}_counts.tsv"
    result_dir="${RESULTS}/${ds}"

    if [ ! -f "${counts_tsv}" ]; then
        echo "  ${ds}: SKIPPING"
        continue
    fi

    echo "  ${ds}: running BayesPrism..."
    Rscript "${SCRIPTS}/12_run_bayesprism.R" \
        "${NEW_SCE}" \
        "${counts_tsv}" \
        "${result_dir}" \
        "${ds}" \
        "celltype" || echo "  WARNING: BayesPrism failed for ${ds}"
done

echo ""
echo "=== Deconvolution Pipeline Complete ==="
echo "End: $(date)"
echo "Results in: ${RESULTS}/"
