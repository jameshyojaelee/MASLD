#!/usr/bin/env Rscript
# KEY MESSAGE: The two co-primary continuum scores identify the same broad, bidirectional transcript-remodeling axis beyond recorded fibrosis stage.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 61_render_transcript_overview.R <analysis_candidate> <figure_candidate>")
}

analysis_candidate <- normalizePath(args[[1]], mustWork = TRUE)
figure_candidate <- args[[2]]
if (dir.exists(figure_candidate)) {
  stop("Refusing to overwrite existing figure candidate: ", figure_candidate)
}

panel_dir <- file.path(figure_candidate, "panels")
source_dir <- file.path(figure_candidate, "source_tables")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

membership <- fread(file.path(analysis_candidate, "bulk", "continuum_membership.tsv.gz"))
stopifnot(nrow(membership) == 23370L)

ordinary <- membership[
  axis_component == FALSE & is.finite(beta_signature_pc1) & is.finite(beta_fixed_projection)
]
ordinary[, display_class := fifelse(
  continuum_associated & beta_fixed_projection > 0,
  "Continuum-associated, positive",
  fifelse(
    continuum_associated & beta_fixed_projection < 0,
    "Continuum-associated, negative",
    "Not continuum-associated"
  )
)]
ordinary[, display_class := factor(
  display_class,
  levels = c(
    "Not continuum-associated",
    "Continuum-associated, negative",
    "Continuum-associated, positive"
  )
)]

fixed_roster <- c("GNMT", "MAT1A", "BICC1")
ordinary[, fixed_display_gene := gene_name %in% fixed_roster]

loo <- fread(file.path(analysis_candidate, "bulk", "signature_loo_gene_meta.tsv"))[
  gene_name %in% c("CYP2C19", "IGFBP7")
]
loo_wide <- dcast(
  loo,
  gene_name + gene_id_base ~ axis_id,
  value.var = "beta"
)
setnames(
  loo_wide,
  c("signature_pc1", "fixed_projection"),
  c("beta_signature_pc1", "beta_fixed_projection")
)
loo_wide[, `:=`(
  display_class = factor(
    "Axis component, leave-one-gene-out",
    levels = "Axis component, leave-one-gene-out"
  ),
  fixed_display_gene = TRUE
)]

plot_data <- ordinary[, .(
  gene_id_base,
  gene_name,
  beta_signature_pc1,
  beta_fixed_projection,
  display_class,
  fixed_display_gene,
  evidence_role = "ordinary_continuum_inference"
)]
fwrite(
  plot_data,
  file.path(source_dir, "continuum_transcript_coprimary_agreement.tsv.gz"),
  sep = "\t"
)
fwrite(
  loo_wide[, .(
    gene_id_base,
    gene_name,
    beta_signature_pc1,
    beta_fixed_projection,
    evidence_role = "axis_component_leave_one_out"
  )],
  file.path(source_dir, "continuum_transcript_axis_components_loo.tsv"),
  sep = "\t"
)

counts <- ordinary[, .N, by = display_class]
counts[, legend_label := paste0(as.character(display_class), " (n=", format(N, big.mark = ","), ")")]
label_map <- setNames(as.character(counts$legend_label), as.character(counts$display_class))

colors <- c(
  "Not continuum-associated" = "#9E9E9E",
  "Continuum-associated, negative" = "#1565C0",
  "Continuum-associated, positive" = "#C9265E"
)

plot <- ggplot(ordinary, aes(beta_signature_pc1, beta_fixed_projection)) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = "#6F6F6F", linewidth = 0.3,
              linetype = "22") +
  geom_point(
    aes(color = display_class),
    size = 0.32,
    alpha = 0.42,
    stroke = 0
  ) +
  geom_point(
    data = ordinary[fixed_display_gene == TRUE],
    shape = 21,
    size = 1.8,
    stroke = 0.35,
    fill = "white",
    color = "black"
  ) +
  geom_point(
    data = loo_wide,
    aes(beta_signature_pc1, beta_fixed_projection),
    inherit.aes = FALSE,
    shape = 23,
    size = 2,
    stroke = 0.4,
    fill = "white",
    color = "black"
  ) +
  geom_text(
    data = ordinary[fixed_display_gene == TRUE],
    aes(label = gene_name),
    size = 6 / ggplot2::.pt,
    family = "Helvetica",
    fontface = "italic",
    nudge_x = 0.028,
    nudge_y = 0.028,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  geom_text(
    data = loo_wide,
    aes(beta_signature_pc1, beta_fixed_projection, label = gene_name),
    inherit.aes = FALSE,
    size = 6 / ggplot2::.pt,
    family = "Helvetica",
    fontface = "italic",
    nudge_x = 0.028,
    nudge_y = 0.028,
    check_overlap = TRUE
  ) +
  scale_color_manual(
    values = colors,
    breaks = names(colors),
    labels = label_map[names(colors)],
    drop = FALSE,
    name = NULL
  ) +
  coord_equal() +
  labs(
    x = "Cohort-PC1 adjusted beta",
    y = "Fixed-projection adjusted beta"
  ) +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    text = element_text(size = 6, face = "plain"),
    axis.text = element_text(size = 6, color = "black", face = "plain"),
    axis.title = element_text(size = 6, face = "plain"),
    legend.text = element_text(size = 6, face = "plain"),
    legend.position = "bottom",
    legend.key.size = grid::unit(3, "mm"),
    legend.box = "vertical",
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1.3, alpha = 1)))

grDevices::pdf.options(useDingbats = FALSE)
ggsave(
  file.path(panel_dir, "s3_continuum_transcript_coprimary_agreement.pdf"),
  plot,
  width = 4.4,
  height = 4.1,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

summary <- data.table(
  metric = c(
    "complete_transcript_family",
    "ordinary_testable_transcripts",
    "continuum_associated_transcripts",
    "axis_components_excluded_from_ordinary_inference",
    "spearman_coprimary_beta_agreement"
  ),
  value = c(
    nrow(membership),
    nrow(ordinary),
    sum(ordinary$continuum_associated),
    sum(membership$axis_component),
    cor(ordinary$beta_signature_pc1, ordinary$beta_fixed_projection, method = "spearman")
  )
)
fwrite(summary, file.path(figure_candidate, "transcript_plot_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(figure_candidate, "sessionInfo.txt"))

