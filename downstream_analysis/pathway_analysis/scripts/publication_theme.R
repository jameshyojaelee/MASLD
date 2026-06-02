# publication_theme.R - Shared color themes for pathway analysis plots
# Defines palettes locally (avoids trailing comma bug in external file)

# Palettes from Sanjana Lab publication standards
# Configure PDF device to use valid fonts (Editable in Illustrator)
pdf.options(useDingbats = FALSE)

palette1 <- c("#e35070", "#e14b9d", "#d358c7", "#b771e6", "#8d8ff4", "#4baeef", "#00c6d9", "#00d4b9", "#30d796", "#86d079", "#bbc26b", "#e1b172")
palette2 <- c("#40b499", "#3fb8b6", "#5db5d0", "#84afe3", "#aca5ed", "#cf9ced", "#ed97e2", "#ff98d0", "#ffa1bb", "#ffb2a6", "#ffc996", "#fee392")

pink_gradient <- c("#f499c1", "#f490bc", "#f388b7", "#f27fb2", "#f176ae",
                   "#f06ca9", "#e35e9a", "#c94d84", "#af3d6e", "#962d59",
                   "#7d1c46", "#650a33")

# Main color: Magenta from palette1
MAIN_COLOR <- "#e14b9d"  # Magenta
MAIN_COLOR_LIGHT <- "#f499c1"  # Light pink
MAIN_COLOR_DARK <- "#7d1c46"  # Dark magenta

# Gradient for continuous scales (low to high significance)
# Using pink gradient for main visualization
SIGNIFICANCE_GRADIENT <- c("#f499c1", "#e14b9d", "#7d1c46")  # Light to dark magenta

# Alternative: Full pink gradient
MAGENTA_GRADIENT <- pink_gradient

# Collection colors using palette1 (magenta-centric)
COLLECTION_COLORS <- c(
    "Hallmark" = palette1[1],   # #e35070 (coral/magenta)
    "GO:BP" = palette1[2],      # #e14b9d (magenta)
    "GO:MF" = palette1[3],      # #d358c7 (purple-magenta)
    "GO:CC" = palette1[4],      # #b771e6 (purple)
    "KEGG" = palette1[5],       # #8d8ff4 (blue-purple)
    "Reactome" = palette1[6]    # #4baeef (blue)
)

# Alternative collection colors using palette2
COLLECTION_COLORS_ALT <- c(
    "Hallmark" = palette2[7],   # #ed97e2 (pink)
    "GO:BP" = palette2[6],      # #cf9ced (purple)
    "GO:MF" = palette2[5],      # #aca5ed (lavender)
    "GO:CC" = palette2[4],      # #84afe3 (light blue)
    "KEGG" = palette2[3],       # #5db5d0 (cyan)
    "Reactome" = palette2[1]    # #40b499 (teal)
)

# Publication-ready theme settings
PUB_THEME <- theme_minimal(base_family = "Helvetica", base_size = 7) +
    theme(
        text = element_text(family = "Helvetica", color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, color = "black"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        panel.grid.minor = element_blank(),
        legend.position = "right"
    )

# Custom dotplot colors for clusterProfiler
# Uses magenta gradient for p-value coloring
get_dotplot_colors <- function() {
    scale_color_gradient(
        low = MAIN_COLOR,      # Bright magenta = significant
        high = MAIN_COLOR_LIGHT,  # Light pink = less significant
        name = "-log10(p.adj)"
    )
}

# Custom fill scale for barplots
get_barplot_fill <- function() {
    scale_fill_gradient(
        low = MAIN_COLOR_LIGHT,
        high = MAIN_COLOR,
        name = "Enrichment Score"
    )
}

# NES color scale (negative to positive)
get_nes_colors <- function() {
    scale_fill_gradient2(
        low = palette1[6],     # Blue for negative NES
        mid = "white",
        high = MAIN_COLOR,     # Magenta for positive NES
        midpoint = 0,
        name = "NES"
    )
}

cat("Publication color theme loaded.\n")
cat(sprintf("Main color (magenta): %s\n", MAIN_COLOR))
