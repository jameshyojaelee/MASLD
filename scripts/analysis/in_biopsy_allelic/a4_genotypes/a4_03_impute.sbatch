#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --array=1-22
# A4 step 3, one autosome per task (PRESPEC_A section 6):
#  panel: 1kGP high-coverage biallelic SNVs with AF >= 0.01 (bref3);
#  standard imputation of the merged RNA calls (Beagle 5.4, local, seed 20260923);
#  leave-gene-out: one Beagle run per gene batch with that batch's exonic RNA sites
#  +-100 kb excluded; output trimmed to the batch genes' +-1 Mb cis windows.
# A gene's cis genotypes for eQTL/allelic analysis come from its batch run.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out a4_03_impute.sbatch <a4_restricted_dir>
set -euo pipefail
out=$1
c=chr$SLURM_ARRAY_TASK_ID
ext=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/allelic_refs
gtf=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/gencode.v49.primary_chr.gtf
src=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/a4_genotypes
mm=/gpfs/commons/home/jameslee/.local/bin/micromamba
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
java_run() { "$mm" run -n allelic_qtl java -Xmx56g "$@"; }
p=$out/panel; i=$out/imputed; b=$out/batches; mkdir -p "$p" "$i" "$b"
map=$ext/beagle/plink_GRCh38_map/chr_in_chrom_field/plink.chr$c.GRCh38.map
gt=$out/merged/$c.rna_calls.vcf.gz

if [ ! -s "$p/$c.common_snv.bref3" ]; then
  bcftools view -v snps -m2 -M2 -i 'INFO/AF >= 0.01' -Oz -o "$p/$c.common_snv.vcf.gz" \
    "$ext/kgp_highcov/1kGP_high_coverage_Illumina.$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  java_run -jar "$ext/beagle/bref3.27Feb25.75f.jar" "$p/$c.common_snv.vcf.gz" > "$p/$c.common_snv.bref3"
  rm "$p/$c.common_snv.vcf.gz"
fi
ref=$p/$c.common_snv.bref3

"$mm" run -n allelic_qtl python "$src/a4_batches.py" "$gtf" "$c" "$b"

if [ ! -s "$i/$c.standard.vcf.gz.tbi" ]; then
  java_run -jar "$ext/beagle/beagle.27Feb25.75f.jar" gt="$gt" ref="$ref" map="$map" out="$i/$c.standard" \
    chrom="$c" nthreads=8 seed=20260923 gp=true
  bcftools index -t "$i/$c.standard.vcf.gz"
fi

for mask in "$b/$c".batch*.mask.tsv; do
  k=$(basename "$mask" .mask.tsv)
  [ -s "$i/$k.cis.vcf.gz.tbi" ] && continue
  bcftools query -R "$mask" -f '%ID\n' "$gt" | sort -u > "$i/$k.excluded_markers.txt"
  java_run -jar "$ext/beagle/beagle.27Feb25.75f.jar" gt="$gt" ref="$ref" map="$map" out="$i/$k.full" \
    chrom="$c" nthreads=8 seed=20260923 gp=true excludemarkers="$i/$k.excluded_markers.txt"
  bcftools index -t "$i/$k.full.vcf.gz"
  bcftools view -R "${mask%.mask.tsv}.cis.tsv" -Oz -o "$i/$k.cis.vcf.gz" "$i/$k.full.vcf.gz"
  bcftools index -t "$i/$k.cis.vcf.gz"
  rm "$i/$k.full.vcf.gz" "$i/$k.full.vcf.gz.tbi"
  echo "$k masked markers $(wc -l < "$i/$k.excluded_markers.txt")"
done
echo "$c done"
