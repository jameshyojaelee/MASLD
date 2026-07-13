BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(BASE, "scripts/figures/publication_theme.R"))
library(ggplot2)
library(dplyr)

# Check if ggrepel is installed, if not we'll use geom_text
use_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (use_ggrepel) {
  library(ggrepel)
}

# Create the dataset based on the user-provided data
drug_data <- data.frame(
  Drug = c("Rezdiffra", "Semaglutide", "FGF21 Analogs", "Survodutide", "Lanifibranor"),
  Stage = c("FDA-approved", "FDA-approved", "Phase 3", "Phase 2b-3", "Phase 3"),
  MASH_Improvement = c(30, 62.9, 30, 55, 49),
  Fibrosis_Improvement = c(25, 36.8, 30, 35, 48),
  Label = c(
    "Rezdiffra (Resmetirom)\nMASH: ~30% [20%]\nFibrosis: 25% [11%]",
    "Semaglutide\nMASH: 62.9% [28.6%]\nFibrosis: 36.8% [14.4%]",
    "FGF21 Analogs\nMASH: ~30% [~28%]\nFibrosis: 30% [~20%]",
    "Survodutide\nMASH: ~55% [~41%]\nFibrosis: ~35% [~13%]",
    "Lanifibranor\nMASH: 49% [27%]\nFibrosis: 48% [19%]"
  )
)

# Convert Stage to a factor to enforce legend order
drug_data$Stage <- factor(drug_data$Stage, levels = c("FDA-approved", "Phase 3", "Phase 2b-3"))

# We will use custom colors for the stages, perhaps drawing from the publication theme
stage_colors <- c(
  "FDA-approved" = "#0D47A1", # Deep blue
  "Phase 3" = "#C2185B",      # Magenta
  "Phase 2b-3" = "#F57F17"    # Warm orange/yellow
)

# Generate the scatter plot
p <- ggplot(drug_data, aes(x = MASH_Improvement, y = Fibrosis_Improvement, color = Stage, shape = Stage)) +
  geom_point(size = 3, stroke = 1) +
  scale_color_manual(values = stage_colors) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(
    x = "MASH Improvement (%)",
    y = "Fibrosis Improvement (%)"
  ) +
  theme_masld() +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "plain")
  ) +
  xlim(15, 75) +
  ylim(15, 55)

# Add text labels
if (use_ggrepel) {
  p <- p + geom_text_repel(aes(label = Label), size = PUB_GEOM_TEXT, show.legend = FALSE, min.segment.length = 0.1, box.padding = 0.5)
} else {
  p <- p + geom_text(aes(label = Label), vjust = -1, size = PUB_GEOM_TEXT, show.legend = FALSE)
}

# Save the plot using the theme's helper function
output_file <- file.path(BASE, "docs/reference/fig_drug_efficacies.pdf")
save_fig(p, output_file, width = fig_full_width, height = 5)

message("[caption] Efficacy of MASLD/MASH Therapeutics: raw efficacy percentages for MASH resolution vs. fibrosis improvement; [#%] = placebo-adjusted efficacy.")
cat("Successfully generated plot at:", output_file, "\n")
