##############################################################################
# Figure 4 (Compact): Pharmacotranscriptomics, Proteomics Validation &
#                      Spatial Context
#
# 8 panels in 3 rows:
#   Row 1: (a) MOA-level LINCS reversal | (b) Clinical drug validation tile
#   Row 2: (c) Ranked enrichment | (d) Effect-size detection | (e) Fibrosis concordance
#   Row 3: (f) Disease-emergent SVGs | (g) Pseudotime dynamics | (h) Enrichment
#
# Sources panel sub-scripts:
#   - fig4_pharma_panels.R     → p_a, p_b
#   - fig4_spatial_panels.R    → p_f, p_g, p_h  (relettered to f,g,h in assembly)
# Proteomics row: 2026-05-28 (P0-H) the former fig4_proteomics_panels.R source
# was removed (deleted 2026-04-23, superseded by fig4_validation.R, run first).
# The canonical proteomics panels are figures/main/fig4_validation/panels/
# fig4{a,b,c}.pdf (the fig4_validation.pdf composite was retired 2026-07-07);
# see the proteomics block below.
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# T1.8 (2026-04-22): fig4_validation.R is the canonical source of reframed
# Fig 4 panels (CHI3L1 SHAP 9.7%, drug/SVG reframed values). Run it in an
# isolated environment FIRST so panels/fig4{a..c}.pdf are rebuilt with the
# canonical values even when the compact wrapper runs. (The fig4_validation.pdf
# composite was retired 2026-07-07 — fig4_validation.R no longer writes it.)
local({
  canonical_script <- file.path(BASE, "scripts/figures/fig4_validation.R")
  if (file.exists(canonical_script)) {
    message("[T1.8] Running canonical fig4_validation.R (reframed panels)...")
    env <- new.env(parent = globalenv())
    env$BASE <- BASE
    tryCatch(
      sys.source(canonical_script, envir = env),
      error = function(e)
        message("WARNING: canonical fig4_validation.R failed: ", e$message)
    )
  } else {
    message("WARNING: canonical fig4_validation.R not found at: ", canonical_script)
  }
})

OUT <- file.path(FIG4_DIR, "fig4_compact.pdf")

# Helper: source a script in isolated environment, extract named panel objects
source_panels <- function(script_path, var_names) {
  panels <- setNames(
    lapply(var_names, function(vn) placeholder(paste0("p_", vn, " [not loaded]"))),
    var_names
  )
  if (!file.exists(script_path)) {
    message("WARNING: Script not found: ", script_path)
    return(panels)
  }
  env <- new.env(parent = globalenv())
  env$BASE <- BASE
  tryCatch({
    source(script_path, local = env)
    for (vn in var_names) {
      var_name <- paste0("p_", vn)
      if (exists(var_name, envir = env))
        panels[[vn]] <- get(var_name, envir = env)
    }
  }, error = function(e) {
    message("WARNING: Error sourcing ", basename(script_path), ": ", e$message)
  })
  panels
}

# ---------------------------------------------------------------------------
# Panels (a-b): Pharmacotranscriptomics
# ---------------------------------------------------------------------------
pharma <- source_panels(
  file.path(BASE, "scripts/figures/fig4_pharma_panels.R"),
  c("a", "b")
)

# ---------------------------------------------------------------------------
# Panels (c-e): Proteomics validation
# 2026-05-28 (P0-H): fig4_proteomics_panels.R was deleted 2026-04-23 and its
# content superseded by the 2026-04-29 rebuild of fig4_validation.R (run first,
# above), which is the canonical proteomics + spatial figure (FIG4_DIR =
# fig4_validation; panels 4a plasma DIA-MS volcano + 4b mRNA-protein
# concordance). The old panels (e/f/g) also relied on GSE276114, which the
# T0.6 correction flagged as RNA-seq mislabeled as proteomics. The source()
# call is removed so the runner no longer silently skips a missing script;
# the canonical proteomics composite is figures/main/fig4_validation/
# fig4_validation.pdf. The compact wrapper marks these slots as a pointer.
proteo <- list(
  e = placeholder("Proteomics: see fig4_validation.pdf (4a/4b)"),
  f = placeholder("Proteomics: see fig4_validation.pdf (4a/4b)"),
  g = placeholder("Proteomics: see fig4_validation.pdf (4a/4b)")
)

# ---------------------------------------------------------------------------
# Panels (f-h): Spatial context
#   Script exports p_f, p_g, p_h — we slot them into figure positions f, g, h
# ---------------------------------------------------------------------------
spatial <- source_panels(
  file.path(BASE, "scripts/figures/fig4_spatial_panels.R"),
  c("f", "g", "h")
)

# ==========================================================================
# Compose: 3 rows — generous sizing to avoid squeezed panels
#   Row 1: 2 pharma panels (wide)
#   Row 2: 3 proteomics panels (proteo e→c, f→d, g→e in figure)
#   Row 3: 3 spatial panels
# ==========================================================================
fig4 <- (pharma[["a"]]  | pharma[["b"]]) /
        (proteo[["e"]]  | proteo[["f"]]  | proteo[["g"]]) /
        (spatial[["f"]] | spatial[["g"]] | spatial[["h"]]) +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 6, face = "plain"))

# Full-width, tall enough that each row has breathing room
save_fig_tall(fig4, OUT, width = fig_full_width, height = 11)
cat("Saved:", OUT, "\n")
