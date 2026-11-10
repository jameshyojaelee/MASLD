#!/bin/bash
#SBATCH --job-name=masld-review-fig7
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/catalog_models/logs/%j.out
#SBATCH --error=Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/catalog_models/logs/%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONNOUSERSITE=1
REVIEW_ROOT=$PWD/Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/catalog_models
RNASEQ_ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
"$RNASEQ_ENV/bin/python" "$REVIEW_ROOT/code/verify_fig7e.py" --project-root "$PWD" --output "$REVIEW_ROOT"
export FIG7_OUT_DIR=$REVIEW_ROOT/fig7/panels FIG7_TABLES=$REVIEW_ROOT/fig7/source_tables FIG7_CODE=$REVIEW_ROOT/fig7/code
"$RNASEQ_ENV/bin/Rscript" --vanilla "$REVIEW_ROOT/fig7/code/fig7e_verdict_vocabulary_and_refusal_thresholds.R"
"$RNASEQ_ENV/bin/Rscript" --vanilla -e 'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("FIG7_OUT_DIR"), "..", "sessionInfo.txt"))'
"$RNASEQ_ENV/bin/python" -m pip freeze > "$REVIEW_ROOT/logs/python_environment.txt"
pdftotext -layout "$REVIEW_ROOT/fig7/panels/fig7E_verdict_vocabulary_and_refusal_thresholds.pdf" "$REVIEW_ROOT/fig7/panels/fig7E.txt"
pdftoppm -scale-to 1800 -png -singlefile "$REVIEW_ROOT/fig7/panels/fig7E_verdict_vocabulary_and_refusal_thresholds.pdf" "$REVIEW_ROOT/fig7/png/fig7E"
