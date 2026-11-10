#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3e/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3e/%x_%j.err
# SS job S3e: composition (InstaPrism, release reference), 7-part ILR, BF4 P1, H38 T2. See PRESPEC_SS_addendum_S3e.md.
# Usage: sbatch s3e_composition.sbatch <S2 dir> <stage tsv> <new out dir> [<ig dir>]
# Mock mode (reproduction): S3E_MOCK=<dataset> builds a release mock into <S2 dir> first (must not exist).
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
IN=$1; STAGE=$2; OUT=$3; IG=${4:-}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
if [[ -n "${S3E_MOCK:-}" ]]; then $RS $SS/scripts/s3e_mock_from_release.R "$S3E_MOCK" "$IN"; fi
$RS $SS/scripts/s3e_composition.R --input "$IN" --stage "$STAGE" --out "$OUT" ${IG:+--ig "$IG"}
$PY -I -m pip freeze > "$OUT/pip_freeze.txt"
sha256sum $SS/scripts/s3e_composition.R $SS/scripts/ss_common.R $SS/scripts/s3e_mock_from_release.R > "$OUT/script_sha256.txt"
echo EXIT_OK
