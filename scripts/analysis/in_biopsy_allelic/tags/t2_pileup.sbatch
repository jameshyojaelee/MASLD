#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=1-22
# T2 input, one autosome per task: allele depths on the 30 Geuvadis BAMs (primary-chr
# alignment, MAPQ 255, BQ 20, depth cap 1,000) at exonic biallelic SNVs polymorphic
# among the same 30 individuals in the 1kGP DNA panel, joined to their DNA genotypes.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%A_%a.out t2_pileup.sbatch <pilot_exec_dir> <out_dir> [excluded_spans.bed]
# The optional BED (t2_exclusion_spans.py: MHC, IG and TR loci) removes those sites before the pileup.
set -euo pipefail
pilot=$1; out=$2; excl=${3:-}
xt=(); [ -n "$excl" ] && xt=(-T "^$excl")
c=chr$SLURM_ARRAY_TASK_ID
w=$pilot/work
ref=/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/genome.primary_chr.fa
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
mkdir -p "$out/t2_work"
bcftools view -i 'AC > 0 && AC < AN' -R "$w/regions/$c.exons.bed" "${xt[@]}" "$w/truth/$c.dna.vcf.gz" -Ou \
  | bcftools query -f '%CHROM\t%POS\t%POS\n' | sort -k2,2n -u > "$out/t2_work/$c.sites.tsv"
bcftools view -i 'AC > 0 && AC < AN' -R "$w/regions/$c.exons.bed" "${xt[@]}" "$w/truth/$c.dna.vcf.gz" -Ou \
  | bcftools query -H -f '%POS\t%REF\t%ALT[\t%GT]\n' > "$out/t2_work/$c.dna.tsv"
bcftools mpileup -f "$ref" -R "$out/t2_work/$c.sites.tsv" -q 255 -Q 20 -d 1000 -a FORMAT/AD \
    -b "$w/calls/bams.txt" -Ou \
  | bcftools query -H -f '%POS\t%REF\t%ALT[\t%AD]\n' > "$out/t2_work/$c.ad.tsv"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python - "$out/t2_work/$c.dna.tsv" "$out/t2_work/$c.ad.tsv" "$out/t2_work/$c.ad_vs_dna.tsv.gz" <<'PY'
import sys, pandas as pd
def clean(cols):
    return [c.split("]", 1)[-1].split(":")[0].split("/")[-1].replace(".Aligned.sortedByCoord.out.bam", "") for c in cols]
dna = pd.read_csv(sys.argv[1], sep="\t"); dna.columns = ["pos", "ref", "alt"] + clean(dna.columns[3:])
ad = pd.read_csv(sys.argv[2], sep="\t"); ad.columns = ["pos", "ref", "alt_pile"] + clean(ad.columns[3:])
rows = []
ad = ad.set_index("pos")
for _, r in dna.iterrows():
    if r["pos"] not in ad.index:
        continue
    a = ad.loc[r["pos"]]
    if isinstance(a, pd.DataFrame):
        a = a.iloc[0]
    alleles = [a["ref"]] + [x for x in str(a["alt_pile"]).split(",") if x != "<*>"]
    if r["alt"] not in alleles:
        continue
    j = alleles.index(r["alt"])
    for s in dna.columns[3:]:
        if s not in a.index:
            continue
        d = [int(v) if v != "." else 0 for v in str(a[s]).split(",")]
        g = str(r[s])
        if "." in g:
            continue
        rows.append((r["pos"], r["ref"], r["alt"], s, d[0], d[j] if j < len(d) else 0, sum(int(q) for q in g.replace("|", "/").split("/"))))
pd.DataFrame(rows, columns=["pos", "ref", "alt", "individual", "ref_reads", "alt_reads", "dna"]).to_csv(sys.argv[3], sep="\t", index=False)
print(len(rows), "site x individual rows")
PY
