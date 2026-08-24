#!/usr/bin/env bash
#SBATCH --job-name=masld-primary-reference
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027

module load SAMtools/1.23.1

package_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
source_fasta="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz"
source_gtf="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
target="${package_root}/executions/reference-build-${SLURM_JOB_ID}"
stage="${target}.staging"
work="$(mktemp -d "${TMPDIR:-/tmp}/masld-primary-reference.XXXXXX")"

cleanup() {
  rm -rf -- "${work}"
  if [[ -d "${stage}" ]]; then
    rm -rf -- "${stage}"
  fi
}
trap cleanup EXIT

if [[ -e "${target}" || -e "${stage}" ]]; then
  printf 'Refusing to overwrite reference output: %s\n' "${target}" >&2
  exit 1
fi
install -d -m 0750 "${stage}"

printf '%s  %s\n' \
  "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215" \
  "${source_fasta}" | sha256sum -c -
printf '%s  %s\n' \
  "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9" \
  "${source_gtf}" | sha256sum -c -
gzip -t "${source_fasta}"
gzip -t "${source_gtf}"

primary_contigs=(
  chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12
  chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
)
printf '%s\n' "${primary_contigs[@]}" > "${work}/primary_contigs.txt"

gzip -dc "${source_fasta}" > "${work}/GRCh38.p14.full.fa"
samtools faidx "${work}/GRCh38.p14.full.fa"
for contig in "${primary_contigs[@]}"; do
  if ! awk -F '\t' -v contig="${contig}" '$1 == contig { found = 1 } END { exit !found }' \
    "${work}/GRCh38.p14.full.fa.fai"; then
    printf 'Required contig is absent from source FASTA: %s\n' "${contig}" >&2
    exit 1
  fi
done

primary_fasta="${stage}/GRCh38.p14.primary_analysis.fa.gz"
samtools faidx \
  --threads 8 \
  --output "${primary_fasta}" \
  --write-index \
  "${work}/GRCh38.p14.full.fa" \
  "${primary_contigs[@]}"
bgzip --test "${primary_fasta}"
test -s "${primary_fasta}.fai"
test -s "${primary_fasta}.gzi"
cut -f 1 "${primary_fasta}.fai" > "${work}/observed_contigs.txt"
cmp "${work}/primary_contigs.txt" "${work}/observed_contigs.txt"
cut -f 1,2 "${primary_fasta}.fai" > "${stage}/GRCh38.p14.primary_analysis.chrom.sizes"

primary_gtf="${stage}/gencode.v49.primary_analysis.annotation.gtf.gz"
awk '
  NR == FNR { allowed[$1] = 1; next }
  /^#/ { print; next }
  $1 in allowed { print }
' "${work}/primary_contigs.txt" <(gzip -dc "${source_gtf}") | bgzip --threads 8 > "${primary_gtf}"
bgzip --test "${primary_gtf}"
gzip -dc "${primary_gtf}" | awk '
  !/^#/ { observed[$1] = 1 }
  END {
    for (contig in observed) {
      if (!(contig ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/)) {
        print "Excluded contig survived GTF filter: " contig > "/dev/stderr"
        exit 1
      }
    }
  }
'

samtools --version > "${stage}/samtools.version.txt"
printf '%s\t%s\n' \
  "source_fasta" "${source_fasta}" \
  "source_fasta_sha256" "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215" \
  "source_gtf" "${source_gtf}" \
  "source_gtf_sha256" "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9" \
  > "${stage}/SOURCE.tsv"
(
  cd "${stage}"
  sha256sum \
    GRCh38.p14.primary_analysis.fa.gz \
    GRCh38.p14.primary_analysis.fa.gz.fai \
    GRCh38.p14.primary_analysis.fa.gz.gzi \
    GRCh38.p14.primary_analysis.chrom.sizes \
    gencode.v49.primary_analysis.annotation.gtf.gz \
    samtools.version.txt \
    SOURCE.tsv \
    > SHA256SUMS
)

chmod 0440 "${stage}"/*
mv "${stage}" "${target}"
printf 'Reference candidate: %s\n' "${target}"
