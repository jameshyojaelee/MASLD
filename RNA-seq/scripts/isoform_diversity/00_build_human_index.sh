#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/build_human_index_%j.out
#SBATCH --error=logs/build_human_index_%j.err

# Phase 1.0 — Build a CLEAN primary-assembly human kallisto index for isoform DTU.
# Subsets the GENCODE v49 comprehensive (chr_patch_hapl_scaff) transcript FASTA to
# the 507,365 primary-assembly transcripts (drops 26,375 scaffold-only alt-haplotype
# near-duplicates that inflate multimapping equivalence classes), then kallisto index.
set -euo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ/RNA-seq/scripts/isoform_diversity"

GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
FULL_FA=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation/RNA-seq/results/kallisto/_index/gencode.v49.transcripts.full.fa.gz

OUT=$PROJ/RNA-seq/results/isoform_diversity/_index
mkdir -p "$OUT"
TAG=v49_primary
PRIM_FA=$OUT/gencode.v49.primary.transcripts.fa.gz
IDX=$OUT/gencode.v49.primary.kidx
WL=$OUT/primary_enst_${TAG}.txt

module purge
module load kallisto/0.51.1

echo "[$(date)] Building tx2gene + primary whitelist from GTF"
python3 make_tx2gene.py --gtf "$GTF" --outdir "$OUT" --tag "$TAG"
echo "[$(date)] whitelist size: $(wc -l < "$WL")"

echo "[$(date)] Subsetting FASTA to primary transcripts (awk; FASTA headers are pipe-delimited: ENST|ENSG|...)"
# Match the versioned ENST (header field 1, before first '|') against the whitelist.
zcat "$FULL_FA" | awk -v wlf="$WL" '
  BEGIN{ while((getline l < wlf) > 0) keep[l]=1 }
  /^>/ { h=substr($0,2); n=split(h,p,"|"); k=(p[1] in keep) }
  { if(k) print }
' | gzip > "$PRIM_FA"
echo "[$(date)] primary FASTA seqs: $(zcat "$PRIM_FA" | grep -c '^>')"

echo "[$(date)] Building kallisto index (k=31)"
kallisto index -i "$IDX" -k 31 "$PRIM_FA"
echo "[$(date)] kallisto inspect:"
kallisto inspect "$IDX" 2>&1 | head -20

echo "[$(date)] DONE. Index: $IDX"
