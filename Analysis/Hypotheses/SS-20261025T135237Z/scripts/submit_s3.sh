#!/bin/bash
# SS: submit S3a-S3h after S2 (run on the login node; sbatch calls only). Writes PRESPEC_SS.sha256 first.
# Order: S3e (composition) -> S3a, S3c, S3d (use the S3e ILR/fractions) -> S3b (needs S3a). S3f/g/h independent.
# All on innovation, job name code-review (user instruction 2026-10-26).
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
cd $SS
[ -s CHECKED_BEFORE_OUTCOME.txt ] || { echo "S2 not finished (no CHECKED_BEFORE_OUTCOME.txt)"; exit 1; }
[ -e PRESPEC_SS.sha256 ] && { echo "PRESPEC_SS.sha256 exists; S3 already submitted?"; exit 1; }
sha256sum PRESPEC_SS.md PRESPEC_SS_addendum_*.md PRESPEC_SS_amendment_*.md > PRESPEC_SS.sha256
chmod a-w PRESPEC_SS.sha256
S2DIR=$SS/out_s2_qc
STAGE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv
TS=$(date -u +%Y%m%dT%H%M%SZ)
S3E=$SS/out_s3e-$TS
SB="env -u SLURM_JOB_ID -u SLURM_JOBID sbatch --parsable --account=innovation --qos=innovation --job-name=code-review --partition=cpu,io"
LOG="--output=$SS/logs/s3/%j.out --error=$SS/logs/s3/%j.err"
je=$($SB $LOG scripts/s3e_composition.sbatch "$S2DIR" "$STAGE" "$S3E" "$SS/out_s1d_ig")
ja=$($SB $LOG --dependency=afterok:$je --export=ALL,S2DIR=$S2DIR,ILR=$S3E/ilr7.tsv scripts/run_s3a.sbatch)
jc=$($SB --dependency=afterok:$je --export=ALL,S2DIR=$S2DIR,STAGE=$STAGE,S3EDIR=$S3E,OUT=$SS/out_s3c-$TS scripts/run_s3c.sbatch)
jd=$($SB --dependency=afterok:$je --export=ALL,S2DIR=$S2DIR,STAGE=$STAGE,S3EDIR=$S3E,OUT=$SS/out_s3d-$TS scripts/run_s3d.sbatch)
jf=$($SB --export=ALL,S2DIR=$S2DIR,STAGE=$STAGE,OUT=$SS/out_s3fgh-$TS scripts/run_s3fgh.sbatch)
printf 'S3e %s %s\nS3a %s (out printed in log)\nS3c %s\nS3d %s\nS3fgh %s\n' $je $S3E $ja $jc $jd $jf | tee logs/s3/submitted_$TS.txt
echo "S3b: submit after S3a with S3ADIR from logs/s3/$ja.out"
