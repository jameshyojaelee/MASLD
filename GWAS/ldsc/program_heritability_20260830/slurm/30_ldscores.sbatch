#!/bin/bash
#SBATCH --job-name=sldsc_l2
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=GWAS/ldsc/program_heritability_20260830/logs/30_ldscores_%j.out
#SBATCH --error=GWAS/ldsc/program_heritability_20260830/logs/30_ldscores_%j.err
#
# Compute LD scores for the 123 annotations against the 1000 Genomes Phase 3
# European panel (489 individuals, 9,997,231 autosomal SNPs), then split the
# result into one file set per annotation.
#
# All 123 annotations go through in a single pass per chromosome. The expensive
# part of the calculation is the pairwise correlation between genotypes inside
# each 1 centimorgan window, and that work is shared across annotations; running
# them separately would repeat it 123 times.
#
# The 1 centimorgan window and the HapMap3 SNP list for --print-snps are the
# settings used to build the baselineLD model itself, so our LD scores and the
# baselineLD LD scores describe the same SNPs on the same scale.
#
# Eight chromosomes run at a time. The largest chromosomes dominate the wall
# clock, so the job finishes in roughly the time chromosome 1 takes on its own
# plus one more round.

set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"

BASE=GWAS/ldsc/program_heritability_20260830
PY=RNA-seq/.mamba/unified_env/bin/python
PLINK=data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink
REF=GWAS/ldsc_ref/sldsc

mkdir -p $BASE/ldscores $BASE/ldscores_split $BASE/logs/l2 $BASE/print_snps

# Build the --print-snps list from baselineLD's OWN LD score SNP column rather
# than from a separately shipped HapMap3 list.
#
# The first attempt used data/gsmap_resource/LDSC_resource/hapmap3_snps/hm.<c>.snp
# and produced 17,261 SNPs on chromosome 22 against baselineLD's 17,489, sharing
# only 15,958. LDSC requires identical SNP columns when concatenating two
# --ref-ld-chr sources, so stage 4 failed on every run with "LD Scores for
# concatenation must have identical SNP columns".
#
# The two lists are both HapMap3 and are not the same HapMap3. Deriving the list
# from the very file the output has to match makes them identical by
# construction instead of by trusting two third-party lists to agree. All 17,489
# baselineLD chromosome 22 SNPs are present in this PLINK panel, so nothing is
# lost by doing it this way.
echo ">> deriving --print-snps from the baselineLD LD score files"
for c in $(seq 1 22); do
  LD=$(find $REF/baselineLD_v2.2 -name "baselineLD.$c.l2.ldscore.gz" -print -quit)
  [[ -n $LD ]] || { echo "baselineLD.$c.l2.ldscore.gz not found" >&2; exit 1; }
  zcat "$LD" | awk 'NR>1{print $2}' > $BASE/print_snps/baselineLD.$c.snp
done
wc -l $BASE/print_snps/baselineLD.22.snp

run_one () {
  local c=$1
  micromamba run -n ldsc ldsc.py \
    --l2 \
    --bfile $PLINK/1000G.EUR.QC.$c \
    --ld-wind-cm 1 \
    --annot $BASE/annot/programs.$c.annot.gz \
    --thin-annot \
    --print-snps $BASE/print_snps/baselineLD.$c.snp \
    --out $BASE/ldscores/programs.$c \
    > $BASE/logs/l2/chr$c.out 2>&1
}
export -f run_one
export BASE PLINK REF

seq 1 22 | xargs -P 8 -I{} bash -c 'run_one {}'

echo ">> checking every chromosome produced LD scores"
for c in $(seq 1 22); do
  [[ -s $BASE/ldscores/programs.$c.l2.ldscore.gz ]] || { echo "chr$c failed"; exit 1; }
done

echo ">> splitting into per-annotation file sets"
seq 1 22 | xargs -P 4 -I{} $PY $BASE/scripts/03_split_ldscores.py {} \
  $BASE/annot $BASE/ldscores $BASE/ldscores_split

# Fail here, in seconds, rather than in stage 4 after another hour. The stage 1
# check asserts that baselineLD's .annot.gz matches the PLINK panel, which is
# true and is a different pair. This asserts the pair that actually has to match
# for --ref-ld-chr concatenation: our LD score SNP column against baselineLD's.
echo ">> asserting our LD score SNP column is identical to baselineLD's"
bad=0
for c in $(seq 1 22); do
  LD=$(find $REF/baselineLD_v2.2 -name "baselineLD.$c.l2.ldscore.gz" -print -quit)
  n_mine=$(zcat $BASE/ldscores/programs.$c.l2.ldscore.gz | awk 'NR>1{print $2}' | sort | md5sum | cut -d' ' -f1)
  n_base=$(zcat "$LD" | awk 'NR>1{print $2}' | sort | md5sum | cut -d' ' -f1)
  if [[ "$n_mine" != "$n_base" ]]; then
    echo "chr$c SNP columns DIFFER between our LD scores and baselineLD" >&2
    bad=1
  fi
done
if [[ $bad -eq 1 ]]; then
  echo "Refusing to hand mismatched LD scores to stage 4." >&2
  exit 1
fi
echo "   all 22 chromosomes: SNP columns identical to baselineLD"

echo "LD SCORES DONE"
