#!/bin/bash
#SBATCH --job-name=H22_h2
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
# H22 job C: stop-rule phase, then all S-LDSC regressions and nulls (05),
# then tau*, paired-jackknife differences, BH families and controls (06).
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_C_h2.sbatch <OUT>
set -euo pipefail
OUT="${1:?out dir required}"
test -s "$OUT/ldsc/packed/l2.22.npy"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$ROOT/scripts/analysis/hypotheses_2026q4/H22_sldsc_trait_difference
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
LENV=/gpfs/commons/home/jameslee/micromamba/envs/ldsc
REF=$ROOT/GWAS/ldsc_ref
cat > "$OUT/ldsc/config.tsv" <<EOF
ldsc_python	$LENV/bin/python
ldsc_py	$LENV/bin/ldsc.py
baseline	$REF/sldsc/baselineLD_v2.2/baselineLD.
weights	$REF/sldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
frq	$REF/sldsc/1000G_Phase3_frq/1000G.EUR.QC.
eur_w_ld	$REF/eur_w_ld_chr/
pdff	$ROOT/GWAS/ldsc/munged/pdff.sumstats.gz
cirr	$ROOT/GWAS/ldsc/munged/ghouse_cirr.sumstats.gz
alt	$OUT/ldsc/sumstats/alt.sumstats.gz
height	$OUT/ldsc/sumstats/height.sumstats.gz
EOF
export OMP_NUM_THREADS=1
N_WORKERS=$(( ${SLURM_CPUS_PER_TASK:-8} - 1 ))  # one ldsc --h2 process per CPU, one CPU spare
$RENV/bin/python -I $HERE/05_run_h2.py "$OUT" "$OUT/ldsc/config.tsv" "$N_WORKERS"
$RENV/bin/python -I $HERE/06_summarize.py "$OUT" "$REF/sldsc/baselineLD_v2.2/baselineLD."
echo "H22 job C done"
