#!/usr/bin/env Rscript

# Regenerate the indexed Figure 4 working-candidate manifest. Figure 4A is
# author-supplied Illustrator artwork and remains absent until provided.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
FIG4 <- file.path(BASE, "figures/main/fig4_singlecell_programs")
panel_paths <- c(
  "panels/fig4b_scrna_umap_embeddable.pdf",
  "panels/fig4c_hotspot_stage_heatmap.pdf",
  "panels/fig4d_communication.pdf",
  "panels/fig4e_bulk_projection_pre_genetics.pdf",
  "panels/fig4f_tf_activity.pdf",
  "panels/supplementary/figs4c_hotspot_stage_heatmap_full.pdf",
  "panels/supplementary/figs4d_communication_heatmap.pdf"
)
source_paths <- list.files(
  file.path(FIG4, "source_tables/current_candidate"), full.names = FALSE
)
relative_paths <- c(panel_paths, file.path("source_tables/current_candidate", source_paths))
paths <- file.path(FIG4, relative_paths)
if (any(!file.exists(paths))) {
  stop("Cannot build Figure 4 manifest; missing: ",
       paste(relative_paths[!file.exists(paths)], collapse = ", "))
}
manifest <- data.table(
  role = c(
    "panel_4b", "panel_4c", "panel_4d", "panel_4e", "panel_4f",
    "supplement_4c_full", "supplement_4d", paste0("source_", source_paths)
  ),
  path = relative_paths,
  size_bytes = as.numeric(file.info(paths)$size),
  sha256 = vapply(paths, digest, character(1), algo = "sha256", file = TRUE,
                  serialize = FALSE),
  state = fifelse(grepl("fig4e", relative_paths), "partial_candidate",
                  "validated_candidate")
)
out <- file.path(FIG4, "manifests/current_candidate_manifest.tsv")
fwrite(manifest, out, sep = "\t", quote = FALSE)
message("[fig4 manifest] wrote ", out, " with ", nrow(manifest), " artifacts")
