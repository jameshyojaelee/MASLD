#!/bin/bash
#SBATCH --job-name=pops
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/67_pops_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/67_pops_%j.err

# PoPS Axis-2 step 2: munge full feature matrix (~57.7k features) + fit PoPS on
# MVP_NAFLD_EUR MAGMA scores + emit which-gene nominations for eQTL-absent loci.
# ADDITIVE. Depends on the ENSGID MAGMA job (66_run_pops_magma.sbatch).
set -eo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd $BASE
mkdir -p GWAS/finemapping/results/seqfunc/logs

POPS=$BASE/tools/pops
# SHARED gpfs staging (was node-local /scratch, invisible from the bigmem compute
# node -> that was the original failure: gzip could not find PoPS.features.txt.gz).
SCR=$BASE/data/external/pops
FEATDIR=$SCR/features_raw
MUNGED=$SCR/features_munged
GENE_ANNOT=$POPS/example/data/utils/gene_annot_jun10.txt
CONTROL=$POPS/example/data/utils/features_jul17_control.txt
MAGMA_PREFIX=$BASE/tools/magma/work/pops/MVP_NAFLD_EUR
OUT_PREDS=$BASE/tools/magma/work/pops/MVP_NAFLD_EUR   # pops writes .preds/.coefs/.marginals
NOM_OUT=$BASE/GWAS/finemapping/results/seqfunc/pops_nominations.tsv
mkdir -p $MUNGED

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

# --- 0. ensure decompressed feature file present -----------------------------
if [ ! -s "$SCR/decomp.done" ] || [ ! -s "$FEATDIR/PoPS.features.txt" ]; then
  echo "[pops] decompressing feature matrix"
  mkdir -p $FEATDIR
  zcat $SCR/PoPS.features.txt.gz > $FEATDIR/PoPS.features.txt
fi
echo "[pops] feature file: $(ls -la $FEATDIR/PoPS.features.txt)"

# --- 1. munge features -------------------------------------------------------
if ! ls $MUNGED/pops_features.mat.0.npy >/dev/null 2>&1; then
  echo "[pops] munging feature matrix (nan_policy=mean)"
  python $POPS/munge_feature_directory.py \
    --gene_annot_path $GENE_ANNOT \
    --feature_dir $FEATDIR \
    --nan_policy mean \
    --save_prefix $MUNGED/pops_features \
    --max_cols 5000
fi
NCHUNK=$(ls $MUNGED/pops_features.mat.*.npy | wc -l)
echo "[pops] munged chunks: $NCHUNK"

# --- 2. run PoPS -------------------------------------------------------------
echo "[pops] fitting PoPS on MVP_NAFLD_EUR MAGMA scores"
python $POPS/pops.py \
  --gene_annot_path $GENE_ANNOT \
  --feature_mat_prefix $MUNGED/pops_features \
  --num_feature_chunks $NCHUNK \
  --magma_prefix $MAGMA_PREFIX \
  --control_features_path $CONTROL \
  --out_prefix $OUT_PREDS \
  --verbose
echo "[pops] preds: $OUT_PREDS.preds"; head -3 $OUT_PREDS.preds; wc -l $OUT_PREDS.preds

# --- 3. which-gene nominations for eQTL-absent loci --------------------------
echo "[pops] building which-gene nominations"
python $BASE/GWAS/finemapping/src/67_pops_whichgene.py \
  --preds $OUT_PREDS.preds \
  --gene_annot $GENE_ANNOT \
  --substrate $BASE/GWAS/finemapping/results/seqfunc/variant_substrate_hg38.tsv \
  --out $NOM_OUT

echo "[pops] done -> $NOM_OUT"
