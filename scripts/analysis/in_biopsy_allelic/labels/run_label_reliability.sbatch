#!/bin/bash
# B-Label: stage-label reliability (labels only). Round 2 (2026-10-15): v1.4 R2 and design inputs, R2 fixed
# (spec v1.3), no restricted-tree read, restricted per-person label tables written to OUT/restricted (chmod 750).
# Submit only through the program helper, for example:
#   P/bin/submit.sh A <run_dir> label "<purpose>" --output=<run_dir>/logs/%j.out run_label_reliability.sbatch \
#     <run_dir> N_BOOT N_GOF FROZEN_CROSSWALK [FROZEN_CROSSWALK_SHA256|-] [RUN_TESTS 0|1]
# FROZEN_CROSSWALK must be a read the lead or user has authorised (RULES.md puts seal-* on the read-denylist).
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=90:00:00
set -euo pipefail
OUT=$1; NBOOT=$2; NGOF=$3; XWALK=$4; XSHA=${5:--}; RUN_TESTS=${6:-0}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
EXE=$ROOT/Analysis/MASLD_Model_Benchmark/executions
R2=$EXE/model-a-gateA-r2-20260930T121832Z          # v1.4 (job 1292)
DESIGN=$EXE/model-a-gateA-design-20260930T121812Z  # v1.4
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
SRC=$ROOT/scripts/analysis/in_biopsy_allelic/labels
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "host $(hostname) start $(date -u +%FT%TZ) job ${SLURM_JOB_ID:-none} cpus ${SLURM_CPUS_PER_TASK:-1}"
mkdir -p "$OUT/restricted"; chmod 750 "$OUT/restricted"
echo "$0 $*" > "$OUT/command_line.txt"
$PY -m pip freeze > "$OUT/pip_freeze.txt"
sha256sum "$SRC/label_reliability.py" "$0" > "$OUT/sha256_script.txt"
if [ "$RUN_TESTS" = "1" ]; then
  cd "$SRC" && $PY -m pytest -q -p no:cacheprovider tests
fi
SHAOPT=()
[ "$XSHA" != "-" ] && SHAOPT=(--frozen-crosswalk-sha256 "$XSHA")

$PY "$SRC/label_reliability.py" \
  --gse193066-placement "$EXE/excluded_cohort_axis_placement-20260829T232000Z/gse193066_donor_keyed_placement.tsv" \
  --gse193066-source "$EXE/kleiner-reliability-20260901T122243Z" \
  --gse193066-series "$ROOT/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/metadata/GSE193066_series_matrix.txt.gz" \
  --gse193066-crosswalk "$ROOT/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-STAGE/PREFLIGHT-v1/audits/gse193066_participant_crosswalk.tsv" \
  --gse193066-n1-pairs "$ROOT/Analysis/Hypotheses/N4-eligibility-corrected-20261014T140730Z/pair_type_table_with_n1.tsv" \
  --pxd-source "$EXE/kleiner-reliability-replication-pxd051911-20260901T211358Z" \
  --pxd-meta "$ROOT/data/PXD051911/meta_data.txt" \
  --frozen-crosswalk "$XWALK" "${SHAOPT[@]}" \
  --kinship-pairs "$EXE/in-biopsy-allelic-a1-full-20260924T025858Z/identity/kinship_related_pairs.tsv" \
  --metadata "$ROOT/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv" \
  --design-individuals "$DESIGN/design_individuals.tsv" \
  --r2-moments "$R2/r2_moments_observed.json" \
  --p2-resamples "$R2/p2_resample_development.tsv.gz" \
  --r2-by-resample "$R2/r2_by_resample.tsv.gz" \
  --codex-lambda "$ROOT/docs/technical/agent_exchange/codex_lambda_20260929/lambda_558.json" \
  --out "$OUT" --restricted-out "$OUT/restricted" --n-boot "$NBOOT" --n-gof "$NGOF" \
  --workers "${SLURM_CPUS_PER_TASK:-1}" && status=0 || status=$?
echo "exit $status end $(date -u +%FT%TZ)"
exit $status
