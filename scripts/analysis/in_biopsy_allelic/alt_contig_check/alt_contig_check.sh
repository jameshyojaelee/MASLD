#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# Read-only check: are genes with copies on GRCh38 ALT/patch contigs undercounted
# in the existing bulk featureCounts (which drop multimappers)?
# For 12 BAMs (3 per cohort), rerun featureCounts with each cohort's original
# flags, once as before and once with -M --fraction. Genes whose name also
# appears on a non-primary contig are "ALT-copy" genes; all other primary genes
# are the comparison. Usage: sbatch --job-name=<n> --output=<out>/logs/%j.out alt_contig_check.sh <out_dir>
set -euo pipefail
out=$1
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
gtf=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
fc=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/featureCounts
py=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
mkdir -p "$out/counts"
"$fc" -v > "$out/environment.txt" 2>&1

declare -A flags=( [GSE130970]="-p -B -s 2" [GSE135251]="-p -B -s 2" [GSE213621]="-p -B -s 2" [GSE193066]="-s 2" )
for c in GSE130970 GSE135251 GSE213621 GSE193066; do
  bams=$(ls "$repo/RNA-seq/Human/Patient_Cohorts/results/$c"/alignments/star/*/*.Aligned.sortedByCoord.out.bam | awk 'NR<=3')
  # shellcheck disable=SC2086
  "$fc" -T 12 ${flags[$c]} -a "$gtf" -o "$out/counts/$c.default.txt" $bams
  # shellcheck disable=SC2086
  "$fc" -T 12 ${flags[$c]} -M --fraction -a "$gtf" -o "$out/counts/$c.multi.txt" $bams
done

"$py" "$repo/scripts/analysis/in_biopsy_allelic/alt_contig_check/summarize.py" \
  --gtf "$gtf" --counts-dir "$out/counts" \
  --degs "$repo/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv" \
  --out "$out"
