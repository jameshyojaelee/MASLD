# Publication theme for MASLD atlas paper
# Nature-style: clean, minimal, consistent
# Magenta/Pink/Blue palette
# Source this in each figure script: source(file.path(BASE, "scripts/figures/publication_theme.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# Illustrator-editable PDFs (no Type3 fonts)
grDevices::pdf.options(useDingbats = FALSE)

# ---------------------------------------------------------------------------
# Shared categorical palettes + soft gradients (palette1/2/3, *_gradient)
# ---------------------------------------------------------------------------
# Loaded once and exposed alongside masld_colors. Palettes 1–3 are the
# softer, equal-luminance categorical sets used as the publication default;
# masld_colors stays available for legacy / semantic mappings.
.pub_palette_path <- file.path(Sys.getenv("HOME"), "publication_color_themes.R")
if (file.exists(.pub_palette_path)) {
  source(.pub_palette_path, local = FALSE)
} else {
  warning("publication_color_themes.R not found at ", .pub_palette_path,
          " — palette1/2/3 will be unavailable.")
}

# ---------------------------------------------------------------------------
# Color palette — Magenta / Pink / Blue
# ---------------------------------------------------------------------------
masld_colors <- list(
  # Primary pair (Liang et al. 2025 palette: deep magenta + warm peach + neutral gray)
  up       = "#C9265E",   # Liang deep magenta (disease/upregulated)
  down     = "#1565C0",   # Deep blue (downregulated; NOT control semantically)
  ns       = "#9E9E9E",   # Liang neutral gray (not significant)

  # Disease states
  control  = "#9E9E9E",   # Liang neutral gray (control — invariant, see feedback_control_gray.md)
  masl     = "#F4A674",   # Liang warm peach (MASL/NAFL — early/mild)
  mash     = "#C9265E",   # Liang deep magenta (MASH/NASH — advanced)
  masld    = "#C9265E",   # Liang deep magenta (MASLD/NAFLD)
  fibrosis = "#A01753",   # Warm deep magenta (Liang-adjacent)

  # Aliases (legacy)
  nafl     = "#F4A674",
  nash     = "#C9265E",
  nafld    = "#C9265E",

  # Sex
  male     = "#1A237E",   # Navy
  female   = "#AD1457",   # Rose magenta

  # Concordance
  conserved = "#00695C",  # Deep teal
  human_enriched = "#7B1FA2",  # Violet
  mouse_specific = "#F4A674",  # Liang warm peach
  discordant     = "#C9265E",  # Liang deep magenta
  not_sig        = "#9E9E9E",  # Liang neutral gray

  # Attribution
  hep_intrinsic  = "#0D47A1",  # Deep blue
  comp_driven    = "#C2185B",  # Magenta
  unmasked       = "#7B1FA2",  # Violet

  # Evidence layers
  deg     = "#1565C0",  # Blue
  mr      = "#C2185B",  # Magenta
  twas    = "#E91E63",  # Pink
  gwas    = "#42A5F5",  # Light blue
  deconv  = "#64B5F6",  # Steel blue
  sex     = "#AD1457",  # Rose
  pathway = "#F48FB1"   # Soft pink
)

# Consensus DEG tier colors
tier_consensus_colors <- c(
  Tier1_HighConfidence = "#880E4F",  # Dark magenta
  Tier2_Moderate       = "#E91E63",  # Bright magenta
  Tier3_Exploratory    = "#42A5F5"   # Light blue
)

# Mouse diet model colors
diet_colors <- c(
  FPC    = "#0D47A1",  # Deep blue
  MCD    = "#C2185B",  # Magenta
  HFD    = "#7B1FA2",  # Violet
  CDAHFD = "#E91E63",  # Pink
  LIDPAD = "#42A5F5"   # Light blue
)

# Attribution class colors
attribution_colors <- c(
  Hepatocyte_intrinsic = "#0D47A1",  # Deep blue
  Composition_driven   = "#C2185B",  # Magenta
  Unmasked             = "#7B1FA2"   # Violet
)

# Sex class colors (interaction-based classification from Script 26 v2)
sex_class_colors <- c(
  Male_biased     = "#1A237E",  # Navy
  Female_biased   = "#AD1457",  # Rose magenta
  Concordant      = "#9E9E9E",  # Grey
  Divergent       = "#E91E63",  # Bright magenta
  # Legacy names (backward compatibility with subsampling/archive figures)
  Male_specific   = "#1A237E",  # Navy (same as Male_biased)
  Female_specific = "#AD1457",  # Rose magenta (same as Female_biased)
  Sex_concordant  = "#00695C",  # Teal (legacy)
  Sex_divergent   = "#E91E63",  # Bright magenta (legacy)
  Not_significant = "#BDBDBD"   # Gray
)

# Concordance category colors (display names)
concordance_colors <- c(
  "Conserved"          = "#00695C",  # Deep teal
  "Human-specific"     = "#7B1FA2",  # Violet
  "Mouse-specific"     = "#F48FB1",  # Soft pink
  "Species-discordant" = "#E91E63",  # Bright magenta
  "Other"              = "#BDBDBD"   # Gray
)

# Multi-evidence tier colors (for legacy compatibility)
tier_colors <- c(
  Tier_A = "#880E4F",
  Tier_B = "#E91E63",
  Tier_C = "#42A5F5"
)

# Fibrosis stage gradient (Fig 2)
fibrosis_stage_colors <- c(
  F0 = "#E3F2FD", F1 = "#90CAF9", F2 = "#42A5F5", F3 = "#1565C0", F4 = "#0D47A1"
)

# Progression colors (Fig 2)
progression_colors <- c(
  Early = "#42A5F5", Late = "#AD1457", Transition = "#7B1FA2"
)

# NAS component colors (Fig 2)
nas_colors <- c(
  Steatosis = "#F57F17", Inflammation = "#C2185B", Ballooning = "#7B1FA2"
)

# Single-cell cell-type palette (17 types, colorblind-friendly; Fig 3)
ct_palette <- c(
  "Hepatocytes"             = "#0D47A1",
  "Cholangiocytes"          = "#1565C0",
  "Endothelial cells"       = "#2E7D32",
  "Fibroblasts"             = "#F57F17",
  "Macrophages"             = "#C2185B",
  "Mono+mono derived cells" = "#E91E63",
  "T cells"                 = "#7B1FA2",
  "B cells"                 = "#9C27B0",
  "Resident NK"             = "#00695C",
  "Circulating NK/NKT"      = "#00897B",
  "Plasma cells"            = "#5D4037",
  "Neutrophils"             = "#FF6F00",
  "cDC1s"                   = "#AD1457",
  "cDC2s"                   = "#D81B60",
  "pDCs"                    = "#6A1B9A",
  "Mig.cDCs"                = "#AB47BC",
  "Basophils"               = "#78909C"
)

# Epigenomic evidence colors (Fig 3 SCENIC+/chromVAR)
epigenomic_colors <- c(
  regulon_up   = "#E91E63",  # Bright magenta (activated in disease)
  regulon_down = "#1565C0",  # Deep blue (suppressed in disease)
  chromvar_sig = "#7B1FA2",  # Violet (significant motif change)
  chromvar_ns  = "#BDBDBD",  # Gray (not significant)
  convergent   = "#00695C"   # Teal (epigenomic-transcriptomic convergent)
)

# Spatial domain colors (Fig 4)
spatial_colors <- c(
  pericentral     = "#C2185B",  # Magenta
  periportal      = "#1565C0",  # Blue
  pan_lobular     = "#7B1FA2",  # Violet
  disease_emergent = "#E91E63", # Pink
  stable_svg      = "#42A5F5",  # Light blue
  disease_lost    = "#BDBDBD"   # Gray
)

# Proteomics validation colors (Fig 4)
proteomics_colors <- c(
  concordant = "#00695C",   # Teal (protein-transcript concordant)
  discordant = "#E91E63",   # Magenta (protein-transcript discordant)
  protein_up = "#C2185B",   # Magenta (protein upregulated)
  protein_dn = "#1565C0"    # Blue (protein downregulated)
)

# Drug evidence strength (Fig 4 validation tile)
drug_evidence_colors <- c(
  Strong   = "#880E4F",  # Dark magenta
  Moderate = "#E91E63",  # Bright magenta
  Absent   = "#E0E0E0"   # Light gray
)

# ---------------------------------------------------------------------------
# Theme
# ---------------------------------------------------------------------------
theme_masld <- function(base_size = 7, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title   = element_text(size = base_size + 1, face = "bold", hjust = 0),
      axis.title   = element_text(size = base_size),
      axis.title.y = element_text(angle = 90, margin = margin(r = 2)),
      axis.title.x = element_text(margin = margin(t = 2)),
      axis.text    = element_text(size = max(base_size - 1, 6), color = "black"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text  = element_text(size = max(base_size - 1, 6)),
      strip.text   = element_text(size = base_size, face = "bold"),
      plot.margin  = margin(3, 3, 3, 3),
      legend.key.size    = unit(0.3, "cm"),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.background   = element_blank(),
      axis.line  = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black")
    )
}

# ---------------------------------------------------------------------------
# Publication font sizes (apply on top of theme_masld via `+ theme_pub()`)
# Equalises text across panels in a multi-panel main figure.
# ---------------------------------------------------------------------------
PUB_AXIS_TEXT  <- 5
PUB_AXIS_TITLE <- 6
PUB_LEGEND     <- 5
PUB_LEGEND_TIT <- 6
PUB_TITLE      <- 7
PUB_SUBTITLE   <- 6
PUB_GEOM_TEXT  <- 1.8     # in-plot text (geom_text) — unitless ggplot "size"
PUB_LEGEND_KEY <- unit(0.18, "cm")

theme_pub <- function() {
  theme(
    plot.title    = element_text(size = PUB_TITLE,    face = "bold"),
    plot.subtitle = element_text(size = PUB_SUBTITLE, color = "gray30"),
    plot.caption  = element_text(size = PUB_SUBTITLE - 1, color = "gray35", hjust = 0),
    axis.title    = element_text(size = PUB_AXIS_TITLE),
    axis.text     = element_text(size = PUB_AXIS_TEXT, color = "black"),
    legend.title  = element_text(size = PUB_LEGEND_TIT, face = "bold"),
    legend.text   = element_text(size = PUB_LEGEND),
    legend.key.size = PUB_LEGEND_KEY,
    strip.text    = element_text(size = PUB_AXIS_TITLE, face = "bold")
  )
}

# ---------------------------------------------------------------------------
# Fig 1 categorical mapping — softer, palette1/3-derived alternatives to the
# deeply saturated masld_colors. Use these for fig1 panels so colours align
# with palette1/2/3; legacy masld_colors stays untouched for fig2+ figures.
# ---------------------------------------------------------------------------
fig1_colors <- list(
  up        = "#C9265E",  # Liang deep magenta (disease/upreg)
  down      = "#518dc9",  # palette3[4]   — medium blue (downreg; NOT control)
  mixed     = "#9b75d6",  # palette3[6]   — violet (multi-set / mixed bin)
  ns        = "#9E9E9E",  # Liang neutral gray
  control   = "#9E9E9E",  # Liang neutral gray (control invariant — feedback_control_gray.md)
  nafl      = "#F4A674",  # Liang warm peach (NAFL/MASL)
  borderline = "#cf9ced", # palette2[6]   — lavender (intermediate)
  nash      = "#C9265E",  # Liang deep magenta
  fibrosis  = "#A01753",  # Liang-adjacent deep magenta
  unknown   = "#9E9E9E"
)

# ---------------------------------------------------------------------------
# Figure dimensions (mm -> inches for ggsave)
# ---------------------------------------------------------------------------
fig_full_width <- 180 / 25.4
fig_half_width <-  88 / 25.4
fig_col_width  <- 120 / 25.4

# ---------------------------------------------------------------------------
# Rasterization helper — keeps PDFs small for large scatterplots
# ---------------------------------------------------------------------------
rasterize_layer <- function(layer, dpi = 300) {
  if (nchar(Sys.getenv("NO_RASTERIZE")) > 0) return(layer)
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::rasterise(layer, dpi = dpi)
  } else {
    layer
  }
}

# ---------------------------------------------------------------------------
# Save helpers
# ---------------------------------------------------------------------------
save_fig <- function(plot, filename, width = fig_full_width, height = 4, dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
  ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, device = pdf_device)
}

save_fig_tall <- function(plot, filename, width = fig_full_width, height = 7, dpi = 300) {
  save_fig(plot, filename, width = width, height = height, dpi = dpi)
}

# Panel label automation (use with patchwork)
auto_tag <- function(p) {
  p + patchwork::plot_annotation(tag_levels = "a")
}

# ---------------------------------------------------------------------------
# Placeholder panel
# ---------------------------------------------------------------------------
placeholder <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 2.5, color = "gray50") +
    theme_void()
}
