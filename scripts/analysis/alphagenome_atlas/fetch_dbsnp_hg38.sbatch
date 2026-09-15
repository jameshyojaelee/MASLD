#!/usr/bin/env bash
#SBATCH --job-name=dbsnp_hg38_fetch
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/dbsnp_fetch_%j.out
set -euo pipefail
D=/gpfs/commons/home/jameslee/reference_genome/dbsnp
cd "$D"
# dbSNP b156+, GRCh38.p14 (GCF_000001405.40). wget only: curl is blocked by cluster egress policy.
wget -c -nv "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
wget -c -nv "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi"
wget -c -nv "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.md5" || true
ls -la "$D"
[ -f GCF_000001405.40.gz.md5 ] && md5sum -c GCF_000001405.40.gz.md5 || echo "NO_MD5_PROVIDED"
echo "DBSNP_FETCH_DONE"
