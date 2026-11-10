#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.err
# SS job S0d: pileup site list (scripts/s0d_targets.py), then a format check: bcftools mpileup -R with the 4-column
# list on one release BAM must reproduce the DIG1 job 62 AD at the 10 job 62 sites. Reads no SteatoSITE data.
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z
OUT=$SS/out_s0d_targets
BCF=/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools
FA=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/genome.primary_chr.fa
sha256sum "$SS/scripts/s0d_targets.py" "$SS/scripts/run_s0d.sbatch" > "$SS/logs/script_sha256_s0d.txt"
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I "$SS/scripts/s0d_targets.py" "$OUT"
# format check on one release sample (GSE130970 first sample of job 62)
read -r s d < <(awk 'NR == 2' "$D/out_rna_genotypes/samples.tsv")
bam=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results/$d/alignments/star/$s/$s.Aligned.sortedByCoord.out.bam
awk -F'\t' '$5 == "dig1_62"' "$OUT/s1b_targets_annotated.tsv" | cut -f1-4 > "$OUT/_check62.tsv"
$BCF mpileup -f "$FA" -R "$OUT/_check62.tsv" -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 -B -d 100000 "$bam##idx##$D/out_rna_genotypes/bai/$s.bai" 2>/dev/null \
  | $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' > "$OUT/_check_new.tsv"
$BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' "$D/out_rna_genotypes/vcf/$s.vcf.gz" > "$OUT/_check_dig1.tsv"
echo "format check sample $s ($d):"
if cmp -s "$OUT/_check_new.tsv" "$OUT/_check_dig1.tsv"; then echo "IDENTICAL to DIG1 job 62"; else echo "DIFFERENT"; diff "$OUT/_check_new.tsv" "$OUT/_check_dig1.tsv" || true; fi
rm "$OUT/_check62.tsv"
cp "$OUT/s1b_targets.tsv" "$SS/scripts/s1b_targets.tsv"
sha256sum "$SS/scripts/s1b_targets.tsv" "$OUT/s1b_targets_annotated.tsv" "$OUT/tnfsf10_ld.tsv" > "$OUT/SHA256SUMS"
wc -l "$SS/scripts/s1b_targets.tsv"
