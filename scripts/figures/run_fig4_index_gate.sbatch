#!/bin/bash
#SBATCH --job-name=fig4_index_gate
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4_index_gate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/fig4_index_gate_%j.err

set -euo pipefail

/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla \
  scripts/figures/update_fig4_manifest.R
python3 scripts/manuscript/validate_resource_scope.py
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla -e '
  library(data.table)
  library(digest)
  root <- "figures/main/fig4_singlecell_programs"
  manifest <- fread(file.path(root, "manifests/current_candidate_manifest.tsv"))
  paths <- file.path(root, manifest$path)
  stopifnot(all(file.exists(paths)))
  observed <- vapply(paths, digest, character(1), algo="sha256", file=TRUE,
                     serialize=FALSE)
  stopifnot(identical(unname(observed), unname(manifest$sha256)))
  cat("FIG4_INDEX_CHECKSUMS\tPASS\t", nrow(manifest), "\n", sep="")

  all_tf <- fread(file.path(root, "source_tables/current_candidate/fig4f_all_tf_activity.tsv"))
  shown <- fread(file.path(root, "source_tables/current_candidate/fig4f_display_tf_activity.tsv"))
  stopifnot(nrow(all_tf) == 1072L, uniqueN(all_tf$tf) == 268L)
  stopifnot(nrow(shown) == 80L, uniqueN(shown$tf) == 20L)
  stopifnot(shown[, all(uniqueN(contrast) == 4L), by=tf]$V1)
  stopifnot(shown[direction == "Increased", all(unique(n_bh) >= 3L), by=tf]$V1)
  stopifnot(shown[direction == "Increased", all(unique(mean_abs_score) >= 5), by=tf]$V1)
  stopifnot(shown[, all(unique(same_direction)), by=tf]$V1)
  stopifnot(setequal(shown[direction == "Reduced", unique(as.character(tf))],
                     c("ZEB2", "HNF4A", "SIX5")))
  cat("FIG4F_SOURCE_CONTRACT\tPASS\t20_TFS\n")
'

pdf="figures/main/fig4_singlecell_programs/panels/fig4f_tf_activity.pdf"
test "$(pdfinfo "${pdf}" | awk '/^Pages:/ {print $2}')" = "1"
if strings "${pdf}" | grep -q '/Subtype /Type3'; then
  echo "Type 3 font detected" >&2
  exit 1
fi
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
