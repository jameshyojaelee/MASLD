#!/usr/bin/env Rscript
# generate_species_contrast_plots.R
# contrasting Human vs Mouse ssGSEA signatures to highlight heterogeneity

# Environment setup
conda_lib <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis/.mamba/pathway_analysis/lib/R/library"
if (dir.exists(conda_lib)) {
    .libPaths(c(conda_lib))
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(grid)
})

# Configuration
BASE_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/pathway_analysis"
OUTPUT_DIR <- file.path(BASE_DIR, "results/refactored")
PLOT_DIR <- file.path(BASE_DIR, "plots/refactored")
THEME_PATH <- file.path(BASE_DIR, "scripts/publication_theme.R")

source(THEME_PATH) # Load PUB_THEME and palette1

# 1. Load Data
cat("Loading ssGSEA results...\n")
ssgsea_df <- readRDS(file.path(OUTPUT_DIR, "all_ssgsea_scores.rds"))

# 2. Define Species Mapping
species_map <- list(
    "Hoang" = "Human",
    "Govaere" = "Human",
    "InHouse_MCD" = "Mouse",
    "Ext_MCD_1" = "Mouse",
    "Ext_MCD_2" = "Mouse"
)

# 3. Calculate Mean Delta Score (Disease - Control) per Species
# 3. Calculate Scores for Contrast
cat("Calculating species-specific scores...\n")


# Human: Calculate Delta (Change) - now that we have labels
human_stats <- ssgsea_df %>%
    filter(Dataset %in% c("Hoang", "Govaere")) %>%
    filter(Condition %in% c("Disease", "Control")) %>%
    group_by(Pathway) %>%
    summarise(
        Mean_Disease = mean(Score[Condition == "Disease"], na.rm=TRUE),
        Mean_Control = mean(Score[Condition == "Control"], na.rm=TRUE),
        Human_Delta = Mean_Disease - Mean_Control,
        .groups = "drop"
    ) %>%
    dplyr::select(Pathway, Human_Delta)


# Mouse: Calculate Delta (Change) - using InHouse labels
mouse_stats <- ssgsea_df %>%
    filter(Dataset == "InHouse_MCD") %>%
    filter(Condition %in% c("Disease", "Control")) %>%
    group_by(Pathway) %>%
    summarise(
        Mean_Disease = mean(Score[Condition == "Disease"], na.rm=TRUE),
        Mean_Control = mean(Score[Condition == "Control"], na.rm=TRUE),
        Mouse_Delta = Mean_Disease - Mean_Control,
        .groups = "drop"
    )

# Join
hallmark_stats <- inner_join(human_stats, mouse_stats, by="Pathway") %>%
    filter(grepl("HALLMARK", Pathway)) %>%
    mutate(
        PathwayName = gsub("HALLMARK_", "", Pathway),
        # Divergence: Human High + Mouse Down (Quadrant IV) vs Human High + Mouse Up (Quadrant I)
        # We define "Diff" technically for ranking, but the plot axes are different units
        # Just use Abs(Mouse_Delta) for coloring or similar
        AbsDelta = abs(Mouse_Delta)
    )

print(head(hallmark_stats))

# 4. Filter for Hallmark
# (Already done in join)

# 5. Plot 1: Species Scatterplot
cat("Generating Scatterplot...\n")

# Highlight specific pathways
highlights <- c("INFLAMMATORY_RESPONSE", "FATTY_ACID_METABOLISM", "TNFA_SIGNALING_VIA_NFKB", "CHOLESTEROL_HOMEOSTASIS")


p_scatter <- ggplot(hallmark_stats, aes(x = Human_Delta, y = Mouse_Delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 1, intercept = 0, color = "grey90", linetype = "dotted") +
    geom_point(aes(color = AbsDelta, size = AbsDelta), alpha = 0.8) +
    geom_text(aes(label = ifelse(PathwayName %in% highlights | AbsDelta > 0.15, PathwayName, "")),
                    size = 3, check_overlap = TRUE, vjust = -1) +
    scale_color_gradient(low = "grey30", high = MAIN_COLOR) +
    annotate("text", x = 0.2, y = -0.15, label = "Discordant\n(Human Up,\n Mouse Down)", color = "blue", alpha = 0.6) +
    annotate("text", x = 0.2, y = 0.15, label = "Conserved\nUpregulated", color = "darkgreen", alpha = 0.6) +
    annotate("text", x = -0.2, y = -0.15, label = "Conserved\nDownregulated", color = "darkgreen", alpha = 0.6) +
    labs(
        title = "Species Contrast: Human NASH vs Mouse MCD (Disease Change)",
        subtitle = "X: Human Delta | Y: Mouse Delta (Disease - Control)",
        x = "Human NASH Delta (Disease - Control)",
        y = "Mouse MCD Delta (Disease - Control)",
        color = "Magnitude"
    ) +
    PUB_THEME +
    theme(legend.position = "none")

ggsave(file.path(PLOT_DIR, "species_contrast_scatterplot.pdf"), p_scatter, width = 8, height = 8)

# 6. Plot 2: Divergence Barplot
cat("Generating Divergence Barplot...\n")

# Divergence Metric: Difference between Deltas
top_divergent <- hallmark_stats %>%
    mutate(Diff = Mouse_Delta - Human_Delta) %>%
    arrange(Diff) %>%
    # Taking top and bottom 10
    slice(c(1:10, (n()-9):n())) %>%
    mutate(Direction = ifelse(Diff > 0, "Mouse Higher / Human Lower", "Mouse Lower / Human Higher"))

p_bar <- ggplot(top_divergent, aes(x = reorder(PathwayName, Diff), y = Diff, fill = Direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Mouse Higher / Human Lower" = "orange", "Mouse Lower / Human Higher" = "purple")) +
    labs(
        title = "Top Divergent Pathways",
        subtitle = "Metric: Mouse Delta - Human Delta",

        y = "Divergence Score",
        x = ""
    ) +
    PUB_THEME

ggsave(file.path(PLOT_DIR, "species_divergence_barplot.pdf"), p_bar, width = 8, height = 6)

cat("Done. Plots saved to", PLOT_DIR, "\n")
