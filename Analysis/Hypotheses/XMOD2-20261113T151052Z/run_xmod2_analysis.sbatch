#!/bin/bash
#SBATCH --job-name=XMOD2an
#SBATCH --account=nslab --qos=nslab --partition=cpu --prefer=v6
#SBATCH --cpus-per-task=4 --mem=32G --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/XMOD2-20261113T151052Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/XMOD2-20261113T151052Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/XMOD2-20261113T151052Z/logs/%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/XMOD2-20261113T151052Z
C5IN=$ROOT/Analysis/Hypotheses/CONF1-20261109T221846Z/out_c5_extract_c1b
GENES=$ROOT/scripts/figures/example_genes.tsv
OUT=$RUN/out_${SLURM_JOB_ID}
mkdir "$OUT"
sha256sum -c "$RUN/PRESPEC.sha256"
sha256sum "$RUN/PRESPEC.md" "$RUN/code/"*.R "$RUN/code/"*.py "$GENES" > "$OUT/sha256_code.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
$ENV/Rscript "$RUN/code/c2_hepatocyte_xmod2.R" "$OUT/trio_c2" trio
$ENV/Rscript "$RUN/code/c5_gse289173_xmod2.R" "$C5IN" "$OUT/trio_c5" trio
$ENV/Rscript "$RUN/code/c2_hepatocyte_xmod2.R" "$OUT/examples_c2" "$GENES"
$ENV/Rscript "$RUN/code/c5_gse289173_xmod2.R" "$C5IN" "$OUT/examples_c5" "$GENES"
PYTHONNOUSERSITE=1 $ENV/python -I "$RUN/code/membership_check.py" "$OUT/fig4d_program_membership.tsv"
