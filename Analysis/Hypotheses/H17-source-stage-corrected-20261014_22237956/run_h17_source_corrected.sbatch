#!/bin/bash
#SBATCH --job-name=H17_stage_corrected
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=logs/H17_stage_corrected_%j.out
#SBATCH --error=logs/H17_stage_corrected_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT="$PWD/Analysis/Hypotheses/H17-source-stage-corrected-20261014_${SLURM_JOB_ID}"
test ! -e "$OUT"; mkdir -p "$OUT/scripts"
HERE="$PWD/scripts/analysis/hypotheses_2026q4/H17_lipid_response_core"
cp "$HERE"/*.R "$OUT/scripts/"
export H17_QUALIFIED_ROSTER="$PWD/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/current_source_diagnosis_donor_roster.tsv"
python3 - "$OUT/scripts/h17_lib.R" <<'PY'
import sys
p=sys.argv[1];s=open(p).read();needle='  donor_meta <- fread(PATHS$donor_meta)'
assert s.count(needle)==1
s=s.replace(needle,needle+'''\n  qualified <- fread(Sys.getenv("H17_QUALIFIED_ROSTER"))
  stopifnot(!anyDuplicated(qualified$donor), all(donor_meta$sample %in% qualified$donor))
  j <- match(donor_meta$sample, qualified$donor)
  donor_meta[, disease_stage_coarse := qualified$disease_stage_coarse[j]]
  donor_meta[, exclude_stage_analysis := qualified$exclude_stage_analysis[j]]''')
open(p,'w').write(s)
PY
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PATH="$ENV/bin:$PATH"
HYP="$PWD/Analysis/Hypotheses/H17-20261009T015231Z"
sha256sum "$HERE"/*.R "$H17_QUALIFIED_ROSTER" > "$OUT/input_sha256_before.tsv"
for STEP in 01_de_and_concordance 02_resampling 03_sensitivity_and_call; do
 "$ENV/bin/Rscript" "$OUT/scripts/$STEP.R" "$HYP" "$OUT"
 echo "Completed $STEP"
 if ls "$OUT"/STOP_* >/dev/null 2>&1; then break; fi
done
sha256sum "$HERE"/*.R "$H17_QUALIFIED_ROSTER" > "$OUT/input_sha256_after.tsv"
cmp "$OUT/input_sha256_before.tsv" "$OUT/input_sha256_after.tsv"
cp "$HERE/run_h17_source_corrected.sbatch" "$OUT/"
