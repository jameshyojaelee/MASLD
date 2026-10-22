#!/bin/bash
#SBATCH --job-name=reads_res
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=logs/reads_extraction/res-%j.out
#SBATCH --error=logs/reads_extraction/res-%j.err

# Download REDIportal v3 (hg38) and build the resource tables for the reads extraction.
# Usage: sbatch build_resources.sbatch RUN_DIR MIN_NSAMPLES TAG [EXISTING_SLIM_TSV_GZ]
set -euo pipefail
RUN_DIR="${1:?run dir required}"
MIN_NS="${2:-50}"
TAG="${3:?tag for the output directory}"
SLIM="${4:-}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
REF=/gpfs/commons/home/jameslee/reference_genome
DL="$RUN_DIR/downloads_untrusted/rediportal_v3"
RES="$RUN_DIR/resources/build_${TAG}"
mkdir -p "$DL"
if [[ ! -s "$DL/TABLE1_hg38_v3.txt.gz" ]]; then
  curl -fsSL --retry 3 -m 3000 -o "$DL/TABLE1_hg38_v3.txt.gz" "http://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v3.txt.gz"
fi
gzip -t "$DL/TABLE1_hg38_v3.txt.gz"
sha256sum "$DL/TABLE1_hg38_v3.txt.gz" | tee "$DL/sha256.txt"
"$PY" -m pip freeze > "$RUN_DIR/resources/pip_freeze_build_${TAG}.txt"
"$PY" -I "$ROOT/scripts/analysis/hypotheses_2026q4/reads_extraction/build_resources.py" \
  "$RES" "$DL/TABLE1_hg38_v3.txt.gz" \
  "$REF/refdata-gex-GRCh38-2024-A/fasta/genome.fa" \
  "$REF/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz" \
  /gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools \
  "$REF/dbsnp/GCF_000001405.40.gz" "$MIN_NS" $SLIM
sha256sum "$RES"/* > "$RES/sha256.txt"
