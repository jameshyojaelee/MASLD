#!/bin/bash
#SBATCH --job-name=sldsc_h2
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=10
#SBATCH --mem=180G
#SBATCH --time=72:00:00
#SBATCH --output=GWAS/ldsc/program_heritability_20260830/logs/40_h2_%j.out
#SBATCH --error=GWAS/ldsc/program_heritability_20260830/logs/40_h2_%j.err
#
# Partitioned heritability, one program at a time.
#
# Each regression asks a single question: after accounting for everything in the
# baselineLD v2.2 model (97 annotations covering coding sequence, conservation,
# regulatory marks, allele frequency, LD structure and background selection),
# and after accounting for the whole pool of genes that were expressed and
# variable enough to enter the Hotspot analysis in that cell type, does this
# particular program's 100 kb windows carry extra heritability?
#
# The background annotation is the reason a program cannot win simply by being
# gene dense: the background covers 50 to 68 percent of the autosome and every
# program's windows sit inside it.
#
# Two more runs per trait provide reference points:
#   _background_only   the background annotation alone on top of baselineLD, so
#                      the background's own enrichment is visible
#   _cts               LDSC's cell-type-specific mode over all 117 programs at
#                      once, a fast independent check on the coefficient
#                      p-values from the per-program runs
#
# Traits marked "excluded" in config/traits.tsv are not run. Their whole-genome
# heritability z-scores or their LD score regression intercepts make partitioned
# estimates uninterpretable; see the README.

set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"

BASE=GWAS/ldsc/program_heritability_20260830
REF=GWAS/ldsc_ref/sldsc
MG=GWAS/ldsc/munged
PY=RNA-seq/.mamba/unified_env/bin/python

# The tarballs decide their own directory layout, so locate one chromosome of
# each and build the prefix from where it landed. An empty result stops the job
# rather than silently pointing LDSC at the working directory.
locate_prefix () {  # locate_prefix <search_root> <example_filename> <prefix_stem>
  local hit
  hit=$(find "$1" -name "$2" | head -1)
  [[ -n $hit ]] || { echo "could not find $2 under $1; run step 1 first" >&2; exit 1; }
  echo "$(dirname "$hit")/$3"
}
BASELINE=$(locate_prefix "$REF/baselineLD_v2.2" 'baselineLD.22.l2.ldscore.gz'  'baselineLD.')
WEIGHTS=$(locate_prefix  "$REF"                 'weights.hm3_noMHC.22.l2.ldscore.gz' 'weights.hm3_noMHC.')
FRQ=$(locate_prefix      "$REF"                 '1000G.EUR.QC.22.frq'          '1000G.EUR.QC.')
SPLIT=$BASE/ldscores_split

echo "baselineLD prefix: $BASELINE"
echo "weights prefix:    $WEIGHTS"
echo "frequency prefix:  $FRQ"

TRAITS=$(awk -F'\t' '!/^#/ && $1!="trait" && NF>1 && $9!="excluded" {print $1}' $BASE/config/traits.tsv)
echo "traits: $TRAITS"

# Build the list of (trait, program) jobs and the cell-type-specific config.
JOBS=$BASE/logs/jobs.txt
: > $JOBS
CTS=$BASE/config/programs.ldcts
: > $CTS
awk -F'\t' 'NR>1 && $3=="program"{print $2"\t"$4}' $BASE/gene_sets/annotation_manifest.tsv \
  | while IFS=$'\t' read -r slug lineage; do
      echo -e "${slug}\t${SPLIT}/${slug}.,${SPLIT}/background-${lineage}." >> $CTS
      for t in $TRAITS; do echo -e "${t}\t${slug}\t${lineage}" >> $JOBS; done
    done
echo "queued $(wc -l < $JOBS) partitioned regressions"

run_one () {
  local trait=$1 slug=$2 lineage=$3
  local out=$BASE/results/$trait
  mkdir -p "$out"
  [[ -s $out/$slug.results ]] && return 0
  micromamba run -n ldsc ldsc.py \
    --h2 $MG/$trait.sumstats.gz \
    --ref-ld-chr "$SPLIT/$slug.,$SPLIT/background-$lineage.,$BASELINE" \
    --w-ld-chr "$WEIGHTS" \
    --frqfile-chr "$FRQ" \
    --overlap-annot \
    --print-coefficients \
    --out "$out/$slug" > /dev/null 2>&1 || echo "FAILED $trait $slug" >&2
}
export -f run_one
export BASE MG SPLIT BASELINE WEIGHTS FRQ

# Reference run: background alone on top of baselineLD.
for t in $TRAITS; do
  mkdir -p $BASE/results/$t
  for lin in hepatocytes fibroblasts macrophages cholangiocytes tcells; do
    micromamba run -n ldsc ldsc.py \
      --h2 $MG/$t.sumstats.gz \
      --ref-ld-chr "$SPLIT/background-$lin.,$BASELINE" \
      --w-ld-chr "$WEIGHTS" --frqfile-chr "$FRQ" \
      --overlap-annot --print-coefficients \
      --out "$BASE/results/$t/_background_only_$lin" > /dev/null 2>&1 \
      || echo "FAILED background $t $lin" >&2
  done
done

echo ">> background annotations on their own, reported before any program result"
$PY $BASE/scripts/06_summarize_backgrounds.py $BASE/results \
    $BASE/results/background_enrichment.tsv

# The 117 program regressions, ten at a time.
awk '{print $1, $2, $3}' $JOBS | xargs -P 10 -n 3 bash -c 'run_one "$0" "$1" "$2"'

# Cell-type-specific mode: an independent check on the coefficient p-values.
for t in $TRAITS; do
  micromamba run -n ldsc ldsc.py \
    --h2-cts $MG/$t.sumstats.gz \
    --ref-ld-chr "$BASELINE" \
    --ref-ld-chr-cts $CTS \
    --w-ld-chr "$WEIGHTS" \
    --out "$BASE/results/$t/_cts" || echo "FAILED cts $t" >&2
done

echo ">> collecting"
$PY $BASE/scripts/05_summarize_results.py $BASE/results $BASE/results/program_heritability_summary.tsv $SPLIT

echo "PARTITIONED HERITABILITY DONE"
