#!/usr/bin/env bash
#SBATCH --job-name=ortho2align
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=ortho2align_work/logs/ortho2align_%j.out
#SBATCH --error=ortho2align_work/logs/ortho2align_%j.err

###############################################################################
# L9 ortho2align — syntenic lncRNA ortholog discovery
#
# Runs ortho2align (Mylarshchikov 2022 BMC Bioinformatics) on all GENCODE v49
# human lncRNAs against mouse (GRCm39 / vM38) via the hg38→mm39 chain file.
# Produces: bestSignificant.annotation.tsv → parsed into L9_ortho2align.tsv
#
# Input BEDs prepared by upstream BED extraction (see task notes).
# Post-processing: parse_ortho2align.py
###############################################################################

set -euo pipefail

# ---- paths ----------------------------------------------------------------
PROJ="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WORK="$PROJ/Cas13_Library_Design/scripts/ortholog_pipeline/ortho2align_work"
OUTDIR="$WORK/full_run"

HUMAN_BED="$WORK/human_lncrna.bed"
MOUSE_GENES_BED="$WORK/mouse_genes.bed"
HUMAN_GENOME="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
MOUSE_GENOME="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCm39-2024-A/fasta/genome.fa"
CHAIN="$PROJ/data/ncrna_conservation/hg38ToMm39.standard_chr.over.chain"

NCORES="${SLURM_CPUS_PER_TASK:-16}"

# ---- env -------------------------------------------------------------------
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate "$PROJ/.envs/ortho2align"

# ---- validate inputs -------------------------------------------------------
for f in "$HUMAN_BED" "$MOUSE_GENES_BED" "$HUMAN_GENOME" "$MOUSE_GENOME" "$CHAIN"; do
    [[ -f "$f" ]] || { echo "MISSING: $f"; exit 1; }
done

echo "=== ortho2align L9 pipeline ==="
echo "Human lncRNAs : $(wc -l < "$HUMAN_BED")"
echo "Mouse genes   : $(wc -l < "$MOUSE_GENES_BED")"
echo "Cores         : $NCORES"
echo "Output dir    : $OUTDIR"
echo "Start time    : $(date)"

# ---- run ortho2align run_pipeline (all-in-one) ----------------------------
mkdir -p "$OUTDIR"

ortho2align run_pipeline \
    -query_genes "$HUMAN_BED" \
    -query_genome "$HUMAN_GENOME" \
    -subject_annotation "$MOUSE_GENES_BED" \
    -subject_genome "$MOUSE_GENOME" \
    -liftover_chains "$CHAIN" \
    -outdir "$OUTDIR" \
    --annotate \
    --fdr \
    -cores "$NCORES" \
    -word_size 6 \
    -min_ratio 0.05 \
    -merge_dist 2000000 \
    -flank_dist 50000 \
    -sample_size 200 \
    -observations 1000 \
    -fitting kde \
    -threshold 0.05 \
    -gapopen 5 \
    -gapextend 2 \
    -timeout 300 \
    -value block_length \
    -function max

echo "=== ortho2align finished ==="
echo "End time: $(date)"

# ---- quick stats -----------------------------------------------------------
if [[ -f "$OUTDIR/stats.txt" ]]; then
    echo ""
    echo "=== stats.txt ==="
    cat "$OUTDIR/stats.txt"
fi

if [[ -f "$OUTDIR/bestSignificant.subject_orthologs.bed" ]]; then
    N_BEST=$(wc -l < "$OUTDIR/bestSignificant.subject_orthologs.bed")
    echo ""
    echo "Best significant orthologs: $N_BEST"
fi

if [[ -f "$OUTDIR/annotation_files/bestSignificant.annotation.tsv" ]]; then
    N_ANN=$(wc -l < "$OUTDIR/annotation_files/bestSignificant.annotation.tsv")
    echo "Annotated orthologs: $((N_ANN - 1))"
fi

echo ""
echo "=== Now run parse_ortho2align.py to produce L9_ortho2align.tsv ==="
