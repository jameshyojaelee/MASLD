#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --array=1-22
# A4 step 2, one autosome per task: merge the cohort calls into one VCF with one
# library per genetic individual (A1 crosswalk; the library with most identity-panel
# calls; possibly mixed libraries dropped), renamed to individual IDs.
# FFPE cohorts (GSE213621, GSE240729): C>T and G>A sites are dropped unless the same
# site is called in a non-FFPE cohort. Sites kept with >= 20 individuals at DP >= 10.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out a4_02_merge.sbatch <a4_restricted_dir> <a1_identity_dir> <a2_dir>
set -euo pipefail
out=$1; idt=$2; a2=$3
c=chr$SLURM_ARRAY_TASK_ID
ffpe=(GSE213621 GSE240729)
nonffpe=(GSE126848 GSE130970 GSE135251 GSE162694 GSE167523 GSE174478 GSE193066 PRJNA512027)
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
m=$out/merged; mkdir -p "$m"

# library choice (same for every chromosome; written once)
if [ ! -s "$m/library_to_individual.tsv" ]; then
  /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python - "$idt/sample_to_individual.tsv" "$a2/excluded_possibly_mixed_libraries.tsv" "$m/library_to_individual.tsv.tmp.$c" <<'PY'
import sys, pandas as pd
x = pd.read_csv(sys.argv[1], sep="\t"); bad = set(pd.read_csv(sys.argv[2], sep="\t")["run"])
x = x[x["individual_id"].notna() & ~x["run"].isin(bad)]
x = x.sort_values(["individual_id", "autosomal_called"], ascending=[True, False]).drop_duplicates("individual_id")
x[["run", "individual_id", "cohort"]].to_csv(sys.argv[3], sep="\t", index=False, header=False)
PY
  mv -n "$m/library_to_individual.tsv.tmp.$c" "$m/library_to_individual.tsv" 2>/dev/null || rm -f "$m/library_to_individual.tsv.tmp.$c"
fi

# FFPE deamination rule
for k in "${nonffpe[@]}"; do
  bcftools query -i '(REF="C" && ALT="T") || (REF="G" && ALT="A")' -f '%CHROM\t%POS\n' "$out/calls/$k.$c.vcf.gz"
done | sort -k2,2n -u > "$m/$c.nonffpe_ctga.tsv"
parts=()
for k in "${nonffpe[@]}"; do parts+=("$out/calls/$k.$c.vcf.gz"); done
for k in "${ffpe[@]}"; do
  bcftools query -i '(REF="C" && ALT="T") || (REF="G" && ALT="A")' -f '%CHROM\t%POS\n' "$out/calls/$k.$c.vcf.gz" \
    | awk 'NR == FNR {seen[$1 ":" $2] = 1; next} !(($1 ":" $2) in seen)' "$m/$c.nonffpe_ctga.tsv" - \
    | sort -k2,2n -u > "$m/$k.$c.drop_ctga.tsv"
  if [ -s "$m/$k.$c.drop_ctga.tsv" ]; then
    bcftools view -T "^$m/$k.$c.drop_ctga.tsv" -Oz -o "$m/$k.$c.ffpe_filtered.vcf.gz" "$out/calls/$k.$c.vcf.gz"
  else
    bcftools view -Oz -o "$m/$k.$c.ffpe_filtered.vcf.gz" "$out/calls/$k.$c.vcf.gz"
  fi
  bcftools index -t "$m/$k.$c.ffpe_filtered.vcf.gz"
  parts+=("$m/$k.$c.ffpe_filtered.vcf.gz")
done

cut -f1 "$m/library_to_individual.tsv" > "$m/$c.keep_runs.txt"
awk '{print $1 "\t" $2}' "$m/library_to_individual.tsv" > "$m/$c.rename.tsv"
bcftools merge -m none --force-samples -Ou "${parts[@]}" \
  | bcftools view -S "$m/$c.keep_runs.txt" --force-samples -Ou \
  | bcftools view -i 'COUNT(FMT/DP >= 10) >= 20' -Oz -o "$m/$c.runs.vcf.gz"
bcftools reheader -s "$m/$c.rename.tsv" -o "$m/$c.rna_calls.vcf.gz" "$m/$c.runs.vcf.gz"
rm "$m/$c.runs.vcf.gz"
bcftools index -t "$m/$c.rna_calls.vcf.gz"
rm "$m"/*."$c".ffpe_filtered.vcf.gz*
echo "$c individuals $(bcftools query -l "$m/$c.rna_calls.vcf.gz" | wc -l) sites $(bcftools view -H "$m/$c.rna_calls.vcf.gz" | wc -l)"
