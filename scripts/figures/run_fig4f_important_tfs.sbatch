#!/bin/bash
#SBATCH --job-name=fig4f_important
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4f_important_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4f_important_%j.err

set -euo pipefail
umask 027

: "${FIGURE_CANDIDATE_ROOT:?Set a fresh candidate output directory}"
: "${STAGE_RELEASE_ROOT:?Set the validated stage-extension directory}"
test -d "${FIGURE_CANDIDATE_ROOT}"
test -d "${STAGE_RELEASE_ROOT}"
test ! -e "${FIGURE_CANDIDATE_ROOT}/figure4/panels/fig4f_tf_activity.pdf"

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"${RNASEQ_ENV}/bin/Rscript" --vanilla \
  "${PROJECT_ROOT}/scripts/figures/tf_convergence_vsF0.R"
"${RNASEQ_ENV}/bin/Rscript" --vanilla -e \
  'writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("FIGURE_CANDIDATE_ROOT"), "sessionInfo.txt"))'

"${RNASEQ_ENV}/bin/Rscript" --vanilla - <<'RS'
library(data.table)
root <- Sys.getenv("FIGURE_CANDIDATE_ROOT")
all_tf <- fread(file.path(root, "source_tables", "fig4f_all_tf_activity.tsv"))
shown <- fread(file.path(root, "source_tables", "fig4f_display_tf_activity.tsv"))
stopifnot(nrow(all_tf) == 1072L, uniqueN(all_tf$tf) == 268L)
stopifnot(nrow(shown) == 80L, uniqueN(shown$tf) == 20L)
stopifnot(shown[, all(uniqueN(contrast) == 4L), by = tf]$V1)
stopifnot(shown[direction == "Increased", all(unique(n_bh) >= 3L), by = tf]$V1)
stopifnot(shown[direction == "Increased", all(unique(mean_abs_score) >= 5), by = tf]$V1)
stopifnot(shown[, all(unique(same_direction)), by = tf]$V1)
stopifnot(setequal(shown[direction == "Reduced", unique(as.character(tf))],
                   c("ZEB2", "HNF4A", "SIX5")))
stopifnot(shown[, all(unique(bh_family_size) == 268L), by = contrast]$V1)
cat("TF source-cardinality checks passed: 268 tested, 20 displayed, 4 contrasts\n")
RS

pdf="${FIGURE_CANDIDATE_ROOT}/figure4/panels/fig4f_tf_activity.pdf"
test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
if strings "${pdf}" | grep -q '/Subtype /Type3'; then
  echo "Type 3 font detected" >&2
  exit 1
fi
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
mkdir -p "${FIGURE_CANDIDATE_ROOT}/proofs"
gs -q -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r200 \
  -sOutputFile="${FIGURE_CANDIDATE_ROOT}/proofs/fig4f_tf_activity.png" "${pdf}"
python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
