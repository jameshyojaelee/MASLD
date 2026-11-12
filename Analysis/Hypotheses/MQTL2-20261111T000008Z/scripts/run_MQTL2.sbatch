#!/bin/bash
#SBATCH --job-name=mqtl2
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/MQTL2-20261111T000008Z/logs/mqtl2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/MQTL2-20261111T000008Z/logs/mqtl2_%j.err
# MQTL2: one job, no array. ALT (Pan-UKBB EUR) x metabolite / lipid-species coloc at 169 ALT loci.
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/MQTL2-20261111T000008Z
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
export PATH=$ENV:$PATH
OUT=$RUN/out_${SLURM_JOB_ID}
test ! -e "$OUT"; mkdir -p "$OUT"
LOCI=$RUN/inputs/alt_loci_169_plus_ctrl.tsv
echo "host $(hostname) | features $(scontrol show node $(hostname -s) | grep -o 'AvailableFeatures=[^ ]*') | $(date)"
sha256sum $RUN/scripts/* $LOCI $RUN/PRESPEC_MQTL2.md \
  $ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z/out_release_vs_disease_loci/loci.tsv \
  $ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2/locus_effects.tsv \
  $ROOT/GWAS/MR_Data/PanUKBB/biomarkers-30620-ALT-both_sexes-irnt.tsv.bgz \
  $ROOT/GWAS/finemapping/data/sumstats/MVP_ALT_EUR_reformatted_hg19.tsv \
  $ROOT/data/external/chen2023_mqtl/chen2023_biomarker_download_list.tsv \
  $ROOT/data/external/chen2023_mqtl/chen2023_accession_trait_map.tsv \
  $ROOT/data/external/lipidqtl/ottensmann2023_accession_trait_map.tsv \
  $RUN/downloads_manifest.tsv \
  $ROOT/data/broadaway_eqtl/hg19ToHg38.over.chain $ROOT/data/broadaway_eqtl/hg38ToHg19.over.chain \
  > $OUT/script_and_input_sha256.txt
# smoke test inside the same job: 1 Chen + 1 Ottensmann file; stop early if it yields no tests
mkdir -p "$OUT/smoke"
MQTL2_TEST_N=1 Rscript $RUN/scripts/02_coloc.R "$OUT/smoke" "$LOCI" 2 > "$OUT/smoke/log.txt" 2>&1 || { cat "$OUT/smoke/log.txt"; exit 1; }
tail -5 "$OUT/smoke/log.txt"
Rscript $RUN/scripts/02_coloc.R "$OUT" "$LOCI" 8
Rscript $RUN/scripts/03_summarise.R "$OUT" "$LOCI"
echo "done $(date)"
