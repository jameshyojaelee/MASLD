#!/bin/bash
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=72:00:00
# Model A reference downloads (user-approved 2026-09-23). wget only (curl HTTPS
# hangs on this cluster). Skips files already complete; writes SHA256 sums.
# Usage: sbatch --job-name=<n> --output=<log>/%j.out download_refs.sbatch
set -euo pipefail
ext=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external
refs=$ext/allelic_refs
kgp=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
mkdir -p "$refs/kgp_highcov" "$refs/rediportal" "$refs/gtex_v8" "$ext/geuvadis/fastq"

# 1) 1kGP high-coverage phased panel, chr1-22 + chrX (v2), with indexes
cd "$refs/kgp_highcov"
for c in $(seq 1 22); do
  f=1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
  wget -c -q "$kgp/$f" "$kgp/$f.tbi"
done
wget -c -q "$kgp/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz" \
           "$kgp/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi"

# 2) REDIportal hg38 v3 editing sites
cd "$refs/rediportal" && wget -c -q http://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v3.txt.gz

# 3) GTEx v8 single-tissue cis-eQTL; keep liver files only
cd "$refs/gtex_v8"
wget -c -q https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar
tar -tf GTEx_Analysis_v8_eQTL.tar | grep -i "Liver" > liver_members.txt
tar -xf GTEx_Analysis_v8_eQTL.tar -T liver_members.txt

# 4) Geuvadis RNA-seq for the DNA-truth gate: 30 unrelated EUR founders chosen by hash
cd "$ext/geuvadis/fastq"
tail -n +2 "$ext/geuvadis/selected_30_eur.tsv" | cut -f5 | tr ';' '\n' | sed 's#^#http://#' > urls.txt
wget -c -q -i urls.txt

for d in "$refs/kgp_highcov" "$refs/rediportal" "$refs/gtex_v8" "$ext/geuvadis/fastq"; do
  (cd "$d" && find . -maxdepth 3 -type f ! -name 'SHA256SUMS*' -newer "$ext/geuvadis/selected_30_eur.tsv" -print0 | xargs -0 sha256sum > "SHA256SUMS_$(date -u +%Y%m%dT%H%M%SZ).txt")
done
echo downloads done
