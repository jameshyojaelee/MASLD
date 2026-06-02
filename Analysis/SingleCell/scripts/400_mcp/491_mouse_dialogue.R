#!/usr/bin/env Rscript
# 491_mouse_dialogue.R
# Independent DIALOGUE run on mouse pseudobulk (requires separate pseudobulk prep on mouse atlas).
# Compare mouse MCP loadings vs human MCP loadings via 1:1 orthologs.

suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
mouse_atlas <- file.path(root, "Analysis/SingleCell/integration/output/mouse/scalesc_mouse_annotated.h5ad")
if (!file.exists(mouse_atlas)) {
  cat(sprintf("[491] mouse atlas not found: %s\n", mouse_atlas))
  quit(status = 0)
}

# Mouse pseudobulk prep not yet built (would need analog of 410/412/413 for mouse).
# Skeleton:
cat("[491] SKELETON: build mouse pseudobulk -> run DIALOGUE -> compare MCPs to human\n")
cat("[491] Next step: write a mouse-specific 412b + run DIALOGUE.run with diet/condition as phenotype\n")
