#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=interactive
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.err
# SS S3a/S3b verification on a mock built from the 774 release (GSE213621 disease units; no SteatoSITE data).
# 0. checkpoint unit test (s3a_test_checkpoint.R); 1. build mock; 2. restart test (8 perms, blocks of 4): uninterrupted run vs run stopped after block 1 and resumed,
# all tables must be identical; a second --resume with every block done must refit nothing and give the same tables; 3. full run, stored norm factors (must reproduce H51 G0 b F3F4 0.7083, F2 0.6148,
# null mean/SD 0.00619/0.0868 and 0.0169/0.121, PC, NC1 slope, per-gene state counts); 4. debug run with RLE recomputed
# on kept genes (H51 G0_rle slopes 0.7088, 0.6160); 5. S3b on the full run.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1  # one BLAS thread per forked worker (as H51)
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
V=$SS/verify_s3a-$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p $V
$R $SS/scripts/s3a_test_checkpoint.R $V/unit_checkpoint
M=$V/mock
$R $SS/scripts/s3a_build_mock.R --out $M
COMMON="--input $M --stage $M/stage_mock.tsv --strata-column acc_q --ilr $M/ilr_mock.tsv --allow-unhashed-prespec"
$R $SS/scripts/s3a_replication.R $COMMON --norm stored --debug --block-size 4 --out $V/ck_uninterrupted
rc=0; $R $SS/scripts/s3a_replication.R $COMMON --norm stored --debug --block-size 4 --stop-after-blocks 1 --out $V/ck_resumed || rc=$?
[ "$rc" -eq 3 ] || { echo "expected exit 3 from the stopped run, got $rc"; exit 1; }
ls $V/ck_resumed/checkpoints
$R $SS/scripts/s3a_replication.R $COMMON --norm stored --debug --block-size 4 --resume --out $V/ck_resumed
$R $SS/scripts/s3a_compare_runs.R $V/ck_uninterrupted $V/ck_resumed $V/restart_test.tsv
$R $SS/scripts/s3a_replication.R $COMMON --norm stored --debug --block-size 4 --resume --out $V/ck_resumed 2>&1 | tee $V/all_done_rerun.log
! grep -q "^.* Wrote " $V/all_done_rerun.log || { echo "all-done rerun refit a block"; exit 1; }
$R $SS/scripts/s3a_compare_runs.R $V/ck_uninterrupted $V/ck_resumed $V/restart_test_all_done.tsv
$R $SS/scripts/s3a_replication.R $COMMON --norm stored --block-size 48 --out $V/s3a_stored
$R $SS/scripts/s3a_replication.R $COMMON --norm rle --debug --out $V/s3a_rle_debug
$R $SS/scripts/s3b_catalog_sidecar.R --s3a $V/s3a_stored --input $M --stage $M/stage_mock.tsv --out $V/s3b --allow-unhashed-prespec
echo "VERIFY_DIR=$V"
