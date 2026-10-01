#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=90:00:00
# STAR 2.7.11b index on the PRIMARY chromosomes only (chr1-22, X, Y, M) of the
# tree-wide GRCh38 FASTA, with GENCODE v49 primary-analysis annotation restricted
# to the same chromosomes. Dropping the 681 ALT/patch/scaffold contigs keeps reads
# in duplicated regions uniquely mapped (see memory: bulk ALT-contig undercount).
# Usage: sbatch --job-name=<n> --output=<log>/%j.out a3_01_star_index.sbatch
set -euo pipefail
ref=/gpfs/commons/home/jameslee/reference_genome
out=$ref/grch38_primary_chr_gencode_v49
env=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
gtf=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/reference-build-21062075/gencode.v49.primary_analysis.annotation.gtf.gz
mkdir -p "$out/star_2.7.11b"
chroms=$(printf 'chr%s ' $(seq 1 22) X Y M)
# shellcheck disable=SC2086
"$env/samtools" faidx "$ref/refdata-gex-GRCh38-2024-A/fasta/genome.fa" $chroms > "$out/genome.primary_chr.fa"
"$env/samtools" faidx "$out/genome.primary_chr.fa"
zcat "$gtf" | awk -F'\t' 'BEGIN {for (i = 1; i <= 22; i++) k["chr" i] = 1; k["chrX"] = 1; k["chrY"] = 1; k["chrM"] = 1}
  /^#/ || ($1 in k)' > "$out/gencode.v49.primary_chr.gtf"
"$env/STAR" --runMode genomeGenerate --runThreadN 16 --genomeDir "$out/star_2.7.11b" \
  --genomeFastaFiles "$out/genome.primary_chr.fa" --sjdbGTFfile "$out/gencode.v49.primary_chr.gtf" \
  --sjdbOverhang 100 --limitGenomeGenerateRAM 80000000000 --outTmpDir "$out/star_tmp"
(cd "$out" && sha256sum genome.primary_chr.fa gencode.v49.primary_chr.gtf star_2.7.11b/SA star_2.7.11b/Genome > SHA256SUMS.txt)
cat > "$out/README.txt" <<EOF
Built $(date -u +%Y-%m-%dT%H:%M:%SZ) by scripts/analysis/in_biopsy_allelic/a3_pilot/a3_01_star_index.sbatch
FASTA: refdata-gex-GRCh38-2024-A/fasta/genome.fa restricted to chr1-22, X, Y, M
GTF: $gtf restricted to the same chromosomes
STAR $("$env/STAR" --version), sjdbOverhang 100
EOF
echo index done
