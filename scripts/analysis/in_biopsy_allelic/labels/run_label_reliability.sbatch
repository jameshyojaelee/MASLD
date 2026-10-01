#!/bin/bash
# B-Label: stage-label reliability (labels only). Job name and log path are given on the sbatch command line.
# Usage: sbatch --job-name=<12xx> --output=<logs>/%j_<12xx>_labels.out run_label_reliability.sbatch \
#          OUT RESTRICTED_OUT N_BOOT N_GOF FROZEN_CROSSWALK RUN_TESTS(0|1)
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=90:00:00
set -euo pipefail
OUT=$1; ROUT=$2; NBOOT=$3; NGOF=$4; XWALK=$5; RUN_TESTS=${6:-0}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
EXE=$ROOT/Analysis/MASLD_Model_Benchmark/executions
RS=$ROOT/RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted
R2=$EXE/model-a-gateA-r2-20260930T120017Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
SRC=$ROOT/scripts/analysis/in_biopsy_allelic/labels
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "host $(hostname) start $(date -u +%FT%TZ) job ${SLURM_JOB_ID:-none} cpus ${SLURM_CPUS_PER_TASK:-1}"
if [ "$RUN_TESTS" = "1" ]; then
  cd "$SRC" && $PY -m pytest -q -p no:cacheprovider tests
fi

$PY "$SRC/label_reliability.py" \
  --gse193066-placement "$EXE/excluded_cohort_axis_placement-20260829T232000Z/gse193066_donor_keyed_placement.tsv" \
  --gse193066-source "$EXE/kleiner-reliability-20260901T122243Z" \
  --gse193066-series "$ROOT/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE193066/metadata/GSE193066_series_matrix.txt.gz" \
  --gse193066-crosswalk "$ROOT/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-STAGE/PREFLIGHT-v1/audits/gse193066_participant_crosswalk.tsv" \
  --pxd-source "$EXE/kleiner-reliability-replication-pxd051911-20260901T211358Z" \
  --pxd-meta "$ROOT/data/PXD051911/meta_data.txt" \
  --frozen-crosswalk "$XWALK" \
  --kinship-pairs "$EXE/in-biopsy-allelic-a1-full-20260924T025858Z/identity/kinship_related_pairs.tsv" \
  --metadata "$ROOT/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv" \
  --design-individuals "$EXE/model-a-gateA-design-20260930T115256Z/design_individuals.tsv" \
  --ancestry "$RS/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv" \
  --r2-moments "$R2/r2_moments_observed.json" \
  --p2-resamples "$R2/p2_resample_development.tsv.gz" \
  --r2-by-resample "$R2/r2_by_resample.tsv.gz" \
  --codex-lambda "$ROOT/docs/technical/agent_exchange/codex_lambda_20260929/lambda_558.json" \
  --out "$OUT" --restricted-out "$ROUT" --n-boot "$NBOOT" --n-gof "$NGOF" \
  --workers "${SLURM_CPUS_PER_TASK:-1}" && status=0 || status=$?
echo "exit $status end $(date -u +%FT%TZ)"
exit $status
