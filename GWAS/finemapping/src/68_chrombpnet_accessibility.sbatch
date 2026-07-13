#!/usr/bin/env bash
#SBATCH --job-name=chrombpnet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chrombpnet_score_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chrombpnet_score_%j.err

# ===========================================================================
# Stage-3 Axis-3 seqfunc: HepG2 ATAC ChromBPNet ACCESSIBILITY MAGNITUDE +
# TF-FOOTPRINT nomination over the eQTL-absent (colocalizes==FALSE) fine-mapped
# credible-set SNVs. Base-resolution, bias-corrected (Kundaje ChromBPNet).
#
# DIRECTION-INDEPENDENT by construction: the Borzoi/AlphaGenome eQTL-direction
# gate FAILED (auROC ~0.53-0.56) and Borzoi ATAC accessibility direction also
# failed (0.537). We report MAGNITUDE (|logFC|, JSD, active-allele-quantile) +
# WHICH-TF-FOOTPRINT (disrupted_tf_motif) + WHICH-GENE. NOT wired into any scored
# convergence channel (circularity w/ COLOC + epigenomic).
#
# Models: ENCODE HepG2 ATAC ChromBPNet, ENCSR291GJU, 5 folds (ENCSR476FKK /
#   ENCFF839MZM), bias-corrected nobias .h5.
# Scorer: kundajelab/variant-scorer (vendored src/chrombpnet_variant_scorer,
#   commit 0e1e341). Env: micromamba chrombpnet_vs (TF 2.15 Keras-2).
# ADDITIVE ONLY (new files); does NOT edit 46d/78/27a or 60-65.
# ===========================================================================
set -uo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"
SRC="${ROOT}/GWAS/finemapping/src"
VS="${SRC}/chrombpnet_variant_scorer"
CBP="${SEQFUNC}/chrombpnet_hepg2"

GENOME=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
CHROMSIZES="${CBP}/hg38.chrom.sizes"
PEAKS="${CBP}/hepg2_atac_peaks.bed"
HITS="${CBP}/motif_hits.bed"
VARLIST="${CBP}/variants_input.chrombpnet.tsv"
GENEMAP="${CBP}/variant_genemap.tsv"
SCOREDIR="${CBP}/scores"
OUT="${SEQFUNC}/chrombpnet_accessibility.tsv"

MAXPEAKS=30000
TOTALSHUF=30000
BATCH=256
SEED=1234

mkdir -p "${SCOREDIR}" "${SEQFUNC}/logs"

echo "=== ChromBPNet HepG2 ATAC accessibility scorer ==="
echo "Node: $(hostname) | Date: $(date)"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || true
echo ""

# Use the env's python by EXPLICIT PATH (micromamba activate does not take effect
# reliably in the batch context -> `python` would fall back to the base env).
ENV=/gpfs/commons/home/jameslee/micromamba/envs/chrombpnet_vs
PY="${ENV}/bin/python"
export PATH="${ENV}/bin:${PATH}"          # so pybedtools finds bedtools/bigBedToBed
export PYTHONUNBUFFERED=1
export PYTHONNOUSERSITE=1
export TF_CPP_MIN_LOG_LEVEL=1
# expose bundled nvidia CUDA/cuDNN wheels to the loader (TF 2.15 [and-cuda])
SITE_PKGS="$("${PY}" -c 'import site;print(site.getsitepackages()[0])')"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done

for f in "${VARLIST}" "${CHROMSIZES}" "${PEAKS}" "${HITS}" "${GENEMAP}"; do
    [ -s "$f" ] || { echo "FATAL: required input missing/empty: $f"; exit 1; }
done

cd "${VS}"

# --- Stage 1: per-fold variant scoring (5 ENCODE folds) --------------------
FOLDS_OK=0
for FOLD in 0 1 2 3 4; do
    MODEL="${CBP}/models/fold_${FOLD}/model.chrombpnet_nobias.fold_${FOLD}.ENCSR291GJU.h5"
    OPREFIX="${SCOREDIR}/fold${FOLD}"
    SCORED="${OPREFIX}.variant_scores.tsv"
    if [ -s "$MODEL" ]; then :; else echo "FATAL: model missing: $MODEL"; exit 1; fi
    if [ -s "$SCORED" ]; then
        echo "[fold ${FOLD}] scores exist, skipping ($SCORED)"; FOLDS_OK=$((FOLDS_OK+1)); continue
    fi
    echo ""; echo "=== [fold ${FOLD}] scoring $(date) ==="
    "${PY}" "${VS}/src/variant_scoring.py" \
        --list        "${VARLIST}" \
        --genome      "${GENOME}" \
        --model       "${MODEL}" \
        --out_prefix  "${OPREFIX}" \
        --chrom_sizes "${CHROMSIZES}" \
        --peaks       "${PEAKS}" \
        --max_peaks   "${MAXPEAKS}" \
        --total_shuf  "${TOTALSHUF}" \
        --batch_size  "${BATCH}" \
        --schema      chrombpnet \
        --random_seed "${SEED}" \
        --no_hdf5
    rc=$?
    if [ $rc -eq 0 ] && [ -s "$SCORED" ]; then FOLDS_OK=$((FOLDS_OK+1)); echo "[fold ${FOLD}] OK"; else
        echo "[fold ${FOLD}] FAILED rc=$rc"; fi
done
echo ""; echo "folds scored OK: ${FOLDS_OK}/5"
[ "${FOLDS_OK}" -ge 1 ] || { echo "FATAL: no fold scored"; exit 1; }

# --- Stage 2: average across folds -----------------------------------------
SCORELIST=""
for FOLD in 0 1 2 3 4; do
    [ -s "${SCOREDIR}/fold${FOLD}.variant_scores.tsv" ] && SCORELIST="${SCORELIST} fold${FOLD}.variant_scores.tsv"
done
echo ""; echo "=== summary across folds ($SCORELIST) $(date) ==="
"${PY}" "${VS}/src/variant_summary_across_folds.py" \
    --score_dir  "${SCOREDIR}" \
    --score_list ${SCORELIST} \
    --out_prefix "${SCOREDIR}/hepg2_atac" \
    --schema     chrombpnet
MEAN="${SCOREDIR}/hepg2_atac.mean.variant_scores.tsv"
[ -s "$MEAN" ] || { echo "FATAL: summary output missing: $MEAN"; exit 1; }

# --- Stage 3: annotate with peak overlap + TF-motif-hit overlap ------------
echo ""; echo "=== annotation (peaks + motif hits) $(date) ==="
"${PY}" "${VS}/src/variant_annotation.py" \
    --list   "${MEAN}" \
    --peaks  "${PEAKS}" \
    --hits   "${HITS}" \
    --schema chrombpnet \
    --out_prefix "${SCOREDIR}/hepg2_atac_annot"
ANNOT="${SCOREDIR}/hepg2_atac_annot.annotations.tsv"
[ -s "$ANNOT" ] || { echo "FATAL: annotation output missing: $ANNOT"; exit 1; }

# --- Stage 4: assemble canonical Axis-3 deliverable ------------------------
echo ""; echo "=== assemble ${OUT} $(date) ==="
"${PY}" "${SRC}/68_chrombpnet_accessibility.py" \
    --annotations "${ANNOT}" \
    --genemap     "${GENEMAP}" \
    --n_folds     "${FOLDS_OK}" \
    --out         "${OUT}"

echo ""; echo "=== DONE $(date) | folds=${FOLDS_OK}/5 | out=${OUT} ==="
