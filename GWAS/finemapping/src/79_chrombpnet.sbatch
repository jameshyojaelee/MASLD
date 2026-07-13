#!/usr/bin/env bash
#SBATCH --job-name=chrombpnet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/79_chrombpnet_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/79_chrombpnet_%j.err

# Phase-5 src/79 caQTL ACCESSIBILITY panel: ChromBPNet HepG2 ATAC signed logfc over
# the Currin caQTL LEADS (direction_labels_caqtl.tsv). Reuses the vendored
# kundajelab/variant-scorer + models/peaks/hits from src/68's chrombpnet_hepg2/.
# Modality-matched (ACCESSIBILITY, NOT gene-expr). ADDITIVE; edits nothing.
set -uo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SEQFUNC="${ROOT}/GWAS/finemapping/results/seqfunc"
SRC="${ROOT}/GWAS/finemapping/src"
VS="${SRC}/chrombpnet_variant_scorer"
CBP="${SEQFUNC}/chrombpnet_hepg2"          # 68's shared assets (models/peaks/hits)
WORK="${SEQFUNC}/chrombpnet_caqtl"         # our isolated working dir
SCOREDIR="${WORK}/scores"

GENOME=/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa
CHROMSIZES="${CBP}/hg38.chrom.sizes"
PEAKS="${CBP}/hepg2_atac_peaks.bed"
HITS="${CBP}/motif_hits.bed"
VARLIST="${WORK}/variants_input.chrombpnet.tsv"
GENEMAP="${WORK}/variant_genemap.tsv"

MAXVARS="${MAXVARS:-4000}"                  # GPU-heavy: cap caQTL leads (report subsample)
MAXPEAKS=30000
TOTALSHUF=30000
BATCH=256
SEED=1234

mkdir -p "${SCOREDIR}" "${SEQFUNC}/logs"
echo "=== src/79 ChromBPNet caQTL accessibility panel ==="; echo "Node: $(hostname) | $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# --- Stage 0: build caQTL VARLIST + genemap (rnaseq env) --------------------
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
set +u; micromamba activate rnaseq; set -u
python "${SRC}/79_score_panels.py" --panel chrombpnet-prep --max-variants "${MAXVARS}"
[ -s "${VARLIST}" ] || { echo "FATAL: VARLIST not built"; exit 1; }

# --- vendored variant-scorer (chrombpnet_vs env) ---------------------------
ENV=/gpfs/commons/home/jameslee/micromamba/envs/chrombpnet_vs
PY="${ENV}/bin/python"; export PATH="${ENV}/bin:${PATH}"
export PYTHONUNBUFFERED=1 PYTHONNOUSERSITE=1 TF_CPP_MIN_LOG_LEVEL=1
SITE_PKGS="$("${PY}" -c 'import site;print(site.getsitepackages()[0])')"
for nvlib in "${SITE_PKGS}"/nvidia/*/lib; do
    [ -d "$nvlib" ] && export LD_LIBRARY_PATH="${nvlib}:${LD_LIBRARY_PATH:-}"
done
cd "${VS}"

FOLDS_OK=0
for FOLD in 0 1 2 3 4; do
    MODEL="${CBP}/models/fold_${FOLD}/model.chrombpnet_nobias.fold_${FOLD}.ENCSR291GJU.h5"
    OPREFIX="${SCOREDIR}/fold${FOLD}"; SCORED="${OPREFIX}.variant_scores.tsv"
    [ -s "$MODEL" ] || { echo "FATAL: model missing $MODEL"; exit 1; }
    if [ -s "$SCORED" ]; then echo "[fold ${FOLD}] exists, skip"; FOLDS_OK=$((FOLDS_OK+1)); continue; fi
    echo "=== [fold ${FOLD}] scoring $(date) ==="
    "${PY}" "${VS}/src/variant_scoring.py" \
        --list "${VARLIST}" --genome "${GENOME}" --model "${MODEL}" \
        --out_prefix "${OPREFIX}" --chrom_sizes "${CHROMSIZES}" --peaks "${PEAKS}" \
        --max_peaks "${MAXPEAKS}" --total_shuf "${TOTALSHUF}" --batch_size "${BATCH}" \
        --schema chrombpnet --random_seed "${SEED}" --no_hdf5
    [ -s "$SCORED" ] && { FOLDS_OK=$((FOLDS_OK+1)); echo "[fold ${FOLD}] OK"; } || echo "[fold ${FOLD}] FAILED"
done
[ "${FOLDS_OK}" -ge 1 ] || { echo "FATAL: no fold scored"; exit 1; }

SCORELIST=""
for FOLD in 0 1 2 3 4; do
    [ -s "${SCOREDIR}/fold${FOLD}.variant_scores.tsv" ] && SCORELIST="${SCORELIST} fold${FOLD}.variant_scores.tsv"
done
"${PY}" "${VS}/src/variant_summary_across_folds.py" \
    --score_dir "${SCOREDIR}" --score_list ${SCORELIST} \
    --out_prefix "${SCOREDIR}/hepg2_atac" --schema chrombpnet
MEAN="${SCOREDIR}/hepg2_atac.mean.variant_scores.tsv"
[ -s "$MEAN" ] || { echo "FATAL: fold-mean missing"; exit 1; }
"${PY}" "${VS}/src/variant_annotation.py" \
    --list "${MEAN}" --peaks "${PEAKS}" --hits "${HITS}" --schema chrombpnet \
    --out_prefix "${SCOREDIR}/hepg2_atac_annot"
ANNOT="${SCOREDIR}/hepg2_atac_annot.annotations.tsv"
[ -s "$ANNOT" ] || { echo "FATAL: annotation missing"; exit 1; }

# --- assemble signed accessibility feature (rnaseq env) --------------------
set +u; micromamba activate rnaseq; set -u
python "${SRC}/79_score_panels.py" --panel chrombpnet-assemble \
    --annotations "${ANNOT}" --n-folds "${FOLDS_OK}"
echo "=== DONE $(date) folds=${FOLDS_OK}/5 ==="
