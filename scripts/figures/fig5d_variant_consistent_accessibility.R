#!/usr/bin/env Rscript

# Figure 5D: variant-consistent open-chromatin context.
# KEY MESSAGE: Colocalization posterior mass overlaps open chromatin across
# several liver lineages, but most gene-study pairs remain only partially
# localized and peak overlap does not assign a target gene.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
grDevices::pdf.options(useDingbats = FALSE)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
CANDIDATE_ROOT <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
if (!nzchar(CANDIDATE_ROOT)) stop("FIGURE_CANDIDATE_ROOT is required")
OUTPUT_DIR <- file.path(CANDIDATE_ROOT, "figure5", "panels")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

ATAC_ROOT <- file.path(
  BASE,
  "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
)
SOURCE <- file.path(
  ATAC_ROOT,
  "genetics/context/genetic_lineage_context_primary_pairs.tsv"
)
FULL_READY <- file.path(ATAC_ROOT, "FULL_READY")
if (!file.exists(SOURCE) || !file.exists(FULL_READY)) {
  stop("Validated ATAC v3 genetic context is incomplete")
}

d <- fread(SOURCE)
expected_lineages <- c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk")
expected_traits <- c("direct_MASLD", "liver_enzyme")
expected_states <- c(
  "replicated_accessible", "source_dependent", "partial", "indeterminate", "untestable"
)
stopifnot(
  nrow(d) == 4080L,
  uniqueN(d[, .(gwas_name, ensembl)]) == 816L,
  uniqueN(d$ensembl) == 462L,
  setequal(d$lineage, expected_lineages),
  setequal(d$trait_class, expected_traits),
  all(d$evidence_state %chin% expected_states),
  all(d$primary_signal_pair == TRUE),
  all(d$mapped_posterior_mass >= 0.95 - 1e-12)
)

state_labels <- c(
  replicated_accessible = "Both sources",
  source_dependent = "One source",
  partial = "Partial mass",
  indeterminate = "No overlap",
  untestable = "Untestable"
)
state_definitions <- c(
  replicated_accessible = "at least 0.5 normalized colocalization posterior mass overlaps peaks in both ATAC sources",
  source_dependent = "at least 0.5 posterior mass overlaps peaks in one ATAC source only",
  partial = "positive peak-overlapping posterior mass below 0.5 in each source",
  indeterminate = "no mapped posterior mass overlaps a peak in either source",
  untestable = "less than 0.95 of posterior mass maps uniquely to hg38"
)
lineage_labels <- c(
  hepatocyte = "Hepatocyte",
  stellate = "Stellate",
  macrophage = "Macrophage",
  cholangiocyte = "Cholangiocyte",
  t_nk = "T/NK"
)
trait_labels <- c(
  direct_MASLD = "MASLD/liver fat",
  liver_enzyme = "Liver enzyme"
)

counts <- d[, .(n_pairs = .N), by = .(trait_class, lineage, evidence_state)]
grid <- CJ(
  trait_class = expected_traits,
  lineage = expected_lineages,
  evidence_state = expected_states,
  unique = TRUE
)
counts <- merge(grid, counts, by = c("trait_class", "lineage", "evidence_state"), all.x = TRUE)
counts[is.na(n_pairs), n_pairs := 0L]
denom <- unique(d[, .(trait_class, gwas_name, ensembl)])[, .(denominator = .N), by = trait_class]
counts <- merge(counts, denom, by = "trait_class", all.x = TRUE)
counts[, `:=`(
  proportion = n_pairs / denominator,
  trait_label = unname(trait_labels[trait_class]),
  lineage_label = unname(lineage_labels[lineage]),
  state_label = unname(state_labels[evidence_state]),
  state_definition = unname(state_definitions[evidence_state]),
  inference_unit = "promoted colocalized gene-study primary signal pair",
  context_metric = "normalized SNP.PP.H4 posterior mass overlapping static open chromatin"
)]
counts[, state_order := match(evidence_state, expected_states)]
setorder(counts, trait_class, lineage, state_order)
counts[, state_order := NULL]
stopifnot(
  nrow(counts) == 50L,
  all(counts[, sum(n_pairs), by = .(trait_class, lineage)]$V1 ==
        counts[, unique(denominator), by = .(trait_class, lineage)]$V1)
)
fwrite(
  counts,
  file.path(DATA_DIR, "fig5d_variant_consistent_accessibility.tsv"),
  sep = "\t", quote = FALSE
)

plot_data <- counts[evidence_state != "untestable"]
plot_data[, lineage_label := factor(
  lineage_label,
  levels = rev(unname(lineage_labels[expected_lineages]))
)]
plot_data[, state_label := factor(
  state_label,
  levels = unname(state_labels[expected_states[1:4]])
)]
plot_data[, facet_label := sprintf(
  "%s (n = %s)", trait_label, format(denominator, big.mark = ",")
)]
facet_levels <- sprintf(
  "%s (n = %s)",
  unname(trait_labels[expected_traits]),
  format(denom[match(expected_traits, trait_class), denominator], big.mark = ",")
)
plot_data[, facet_label := factor(facet_label, levels = facet_levels)]

fills <- c(
  "Both sources" = "#00695C",
  "One source" = "#7B1FA2",
  "Partial mass" = "#8EC7E2",
  "No overlap" = "#D9D9D9"
)

p <- ggplot(plot_data, aes(x = proportion, y = lineage_label, fill = state_label)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.12) +
  geom_text(
    aes(label = ifelse(proportion >= 0.065 & n_pairs > 0L, n_pairs, "")),
    position = position_stack(vjust = 0.5),
    family = "Helvetica", size = 6 / .pt, colour = "#202124"
  ) +
  facet_wrap(~facet_label, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = fills, drop = FALSE, name = NULL) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1), labels = c("0", "50", "100"),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(x = "Gene–study pairs (%)", y = NULL) +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    text = element_text(size = 6, face = "plain"),
    axis.text = element_text(size = 6, colour = "black", face = "plain"),
    axis.title = element_text(size = 6, face = "plain"),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    strip.text = element_text(size = 6, face = "plain", hjust = 0),
    strip.background = element_blank(),
    panel.spacing.x = grid::unit(7, "pt"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 6, face = "plain"),
    legend.key.size = grid::unit(0.18, "cm"),
    legend.spacing.x = grid::unit(2, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(1, 5, 1, 1)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

out <- file.path(OUTPUT_DIR, "fig5d_snatac_accessibility.pdf")
ggsave(out, p, width = 2.96, height = 1.88, device = grDevices::cairo_pdf)
message("[fig5d] saved variant-consistent accessibility context: ", out)
message(
  "CAPTION (Fig. 5D): Normalized SNP.PP.H4 posterior mass for the primary signal pair of each ",
  "promoted colocalized gene-study pair was intersected with static open chromatin in two liver ",
  "snATAC sources across five shared lineages. Counts show posterior-mass evidence states, with ",
  "direct MASLD/liver-fat and liver-enzyme traits kept separate. Peak overlap localizes regulatory ",
  "DNA but does not assign a target gene; partial or absent overlap is unresolved context, not a ",
  "negative functional test."
)
