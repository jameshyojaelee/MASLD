
# theme_publication.R
# Sanjana Lab Publication Theme Standard
# Updated: Feb 2026

library(ggplot2)
library(grid)

# 1. Colors
# ---------------------------------------------------------------------------
sanjana_colors <- c(
  Magenta  = "#e14b9d", # Primary Highlight
  Pink     = "#e35070", # Secondary
  Purple   = "#d358c7", # Tertiary
  Blue     = "#4baeef", # Background
  Orange   = "#e1b172", # Contrast/Warning
  Green    = "#30d796", # Contrast/Positive
  Grey     = "#b0b0b0", # Non-significant (Light)
  DarkGrey = "#606060"  # Non-significant (Dark)
)

# Standard condition colors
condition_colors <- c(
  "Control"       = "#4baeef",       # Blue
  "Control_Obese" = "#89CFF0",       # Lighter blue
  "NAFL"          = "#e1b172",       # Orange
  "NASH"          = "#e14b9d",       # Magenta
  "NASH_Fibrosis" = "#d358c7",       # Purple
  "Fibrosis"      = "#d358c7",
  "Fibrosis_F0F1" = "#e1b172",
  "Fibrosis_F1"   = "#e1b172",
  "Fibrosis_F2"   = "#e35070",       # Pink
  "Fibrosis_F2F3" = "#e14b9d",
  "Fibrosis_F3F4" = "#d358c7",
  "Fibrosis_F4"   = "#d358c7"
)

# 2. ggplot2 Theme
# ---------------------------------------------------------------------------
# "Helvetica" should be installed. If not, fallback to "sans".
font_family <- "Helvetica"

theme_publication <- function(base_size = 7, base_family = font_family) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Text
      plot.title    = element_text(size = 8, face = "bold", hjust = 0),
      axis.title    = element_text(size = 8),
      axis.text     = element_text(size = 6),
      legend.text   = element_text(size = 6),
      legend.title  = element_text(size = 7, face = "bold"),
      strip.text    = element_text(size = 7, face = "bold"),
      
      # Grid & Background
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      
      # Axes
      axis.line = element_line(color = "black", linewidth = 0.3),
      
      # Legend
      legend.key.size = unit(0.8, "lines"),
      legend.key = element_rect(fill = NA, color = NA)
    )
}

# 3. PDF Device Helper
# ---------------------------------------------------------------------------
# Wrapper to save PDF with correct settings (no Dingbats, cairo if available)
save_pdf <- function(plot, filename, width = 4, height = 3) {
  # Use cairo_pdf if available for better font handling, else standard pdf
  if (capabilities("cairo")) {
    ggsave(filename, plot = plot, device = cairo_pdf, width = width, height = height)
  } else {
    ggsave(filename, plot = plot, device = "pdf", width = width, height = height, useDingbats = FALSE)
  }
  # No PNGs per guidelines
}
