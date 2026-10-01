#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=90:00:00
#SBATCH --array=1-22
# A3 pilot imputation, one autosome per task. Beagle 5.4 (27Feb25.75f), 1kGP
# high-coverage panel without the 30 test individuals' families, GRCh38 PLINK map.
#   standard: all RNA-called sites as input;
#   masked:   exonic RNA sites of a hash-chosen 10% of genes, +-100 kb, excluded
#             (the leave-gene-out rule of PRESPEC_A), so imputation at those sites
#             uses only sites outside the gene window.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out a3_04_impute.sbatch <exec_dir>
set -euo pipefail
exec_dir=$1
c=chr$SLURM_ARRAY_TASK_ID
w=$exec_dir/work; mkdir -p "$w/imputed"
ext=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/allelic_refs
refdir=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49
mm=/gpfs/commons/home/jameslee/.local/bin/micromamba
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
beagle() { "$mm" run -n allelic_qtl java -Xmx34g -jar "$ext/beagle/beagle.27Feb25.75f.jar" "$@"; }
map=$ext/beagle/plink_GRCh38_map/chr_in_chrom_field/plink.chr$c.GRCh38.map

# masked gene windows: genes on this chromosome ranked by sha256("MASLD-A-mask" + gene_id); lowest 10%
awk -F'\t' -v c="$c" '$1 == c && $3 == "gene" {match($9, /gene_id "[^"]+"/); g = substr($9, RSTART + 9, RLENGTH - 10); print g "\t" $4 "\t" $5}' \
  "$refdir/gencode.v49.primary_chr.gtf" > "$w/imputed/$c.genes.tsv"
"$mm" run -n allelic_qtl python - "$w/imputed/$c.genes.tsv" "$w/imputed/$c.mask_windows.tsv" <<'EOF'
import hashlib, sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
rows.sort(key=lambda r: hashlib.sha256(("MASLD-A-mask" + r[0]).encode()).hexdigest())
keep = rows[: max(1, len(rows) // 10)]
with open(sys.argv[2], "w") as fh:
    for g, s, e in keep:
        fh.write(f"{g}\t{max(1, int(s) - 100000)}\t{int(e) + 100000}\n")
EOF
awk -v c="$c" 'BEGIN {OFS = "\t"} {print c, $2, $3}' "$w/imputed/$c.mask_windows.tsv" | sort -k2,2n > "$w/imputed/$c.mask_regions.tsv"
bcftools query -R "$w/imputed/$c.mask_regions.tsv" -f '%ID\n' "$w/calls/$c.rna.vcf.gz" | sort -u > "$w/imputed/$c.excluded_markers.txt"

beagle gt="$w/calls/$c.rna.vcf.gz" ref="$w/panel/$c.ref.vcf.gz" map="$map" out="$w/imputed/$c.standard" \
  chrom="$c" nthreads=4 seed=20260923 gp=true
beagle gt="$w/calls/$c.rna.vcf.gz" ref="$w/panel/$c.ref.vcf.gz" map="$map" out="$w/imputed/$c.masked" \
  chrom="$c" nthreads=4 seed=20260923 gp=true excludemarkers="$w/imputed/$c.excluded_markers.txt"
bcftools index -t "$w/imputed/$c.standard.vcf.gz"
bcftools index -t "$w/imputed/$c.masked.vcf.gz"
echo "$c masked markers $(wc -l < "$w/imputed/$c.excluded_markers.txt")"
