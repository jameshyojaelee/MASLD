#!/bin/bash
#SBATCH --job-name=A1_h2
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=72:00:00
# A1 job C (adapted from H22 job C): testability phase, then every S-LDSC
# regression including the seeded Q2 nulls (05), then tau*, differences, BH
# families, null calibration and controls (06).
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_C_h2.sbatch <OUT>
set -euo pipefail
OUT="${1:?out dir required}"
test -s "$OUT/ldsc/packed/l2.1.npy"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
HERE=$OUT/scripts
RENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
LENV=/gpfs/commons/home/jameslee/micromamba/envs/ldsc
REF=$ROOT/GWAS/ldsc_ref
H22=$ROOT/Analysis/Hypotheses/H22-20261009T015913Z
cat > "$OUT/ldsc/config.tsv" <<EOF
ldsc_python	$LENV/bin/python
ldsc_py	$LENV/bin/ldsc.py
baseline	$REF/sldsc/baselineLD_v2.2/baselineLD.
weights	$REF/sldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
frq	$REF/sldsc/1000G_Phase3_frq/1000G.EUR.QC.
h22	$H22
EOF
export OMP_NUM_THREADS=1
N_WORKERS=$(( ${SLURM_CPUS_PER_TASK:-16} - 1 ))
$RENV/bin/python -I $HERE/05_run_h2.py "$OUT" "$OUT/ldsc/config.tsv" "$OUT/traits.tsv" "$N_WORKERS"
$RENV/bin/python -I $HERE/06_summarize.py "$OUT" "$REF/sldsc/baselineLD_v2.2/baselineLD." "$OUT/traits.tsv" "$H22"
echo "A1 job C done"
