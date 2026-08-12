#!/usr/bin/env Rscript

# KEY MESSAGE: Colocalized variant posterior mass is observable in zero, one,
# or both ATAC cohorts across multiple liver lineages without forcing one cell type.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(ggplot2)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
RELEASE_ID <- "atac-context-v3-candidate-2026-08-11-r1"
CANDIDATE <- Sys.getenv(
  "ATAC_V3_CANDIDATE_ROOT",
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID)
)
EXPECTED <- normalizePath(
  file.path(BASE, "Analysis/Multimodal_Program_Projection/candidates", RELEASE_ID),
  mustWork = FALSE
)
if (!identical(normalizePath(CANDIDATE, mustWork = FALSE), EXPECTED)) stop("Unsafe candidate root")
SOURCE <- file.path(CANDIDATE, "genetics/context/genetic_lineage_context_primary_pairs.tsv")
if (!file.exists(SOURCE)) stop("Primary variant-consistent lineage table is absent")
OUT <- file.path(CANDIDATE, "genetics/figures")
if (dir.exists(OUT)) stop("Refusing to overwrite genetic figure output")
dir.create(file.path(OUT, "source_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "panels"), recursive = TRUE, showWarnings = FALSE)
grDevices::pdf.options(useDingbats = FALSE)
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))

context <- fread(SOURCE)
counts <- context[, .(n_pairs = .N), by = .(trait_class, lineage, evidence_state)]
denominators <- unique(context[, .(trait_class, lineage, gwas_name, ensembl)])[,
  .(denominator = .N), by = .(trait_class, lineage)
]
counts <- merge(counts, denominators, by = c("trait_class", "lineage"), all.x = TRUE)
setorder(counts, trait_class, lineage, evidence_state)
source_path <- file.path(OUT, "source_tables", "fig4d_variant_consistent_lineage_counts.tsv")
fwrite(counts, source_path, sep = "\t")

counts[, lineage := factor(
  lineage,
  levels = c("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"),
  labels = c("Hepatocyte", "Stellate", "Macrophage", "Cholangiocyte", "T/NK")
)]
counts[, evidence_state := factor(
  evidence_state,
  levels = c("replicated_accessible", "source_dependent", "partial", "indeterminate", "untestable")
)]
plot <- ggplot(counts, aes(lineage, evidence_state, size = n_pairs, color = evidence_state)) +
  geom_point(alpha = 0.85) +
  geom_text(aes(label = paste0(n_pairs, "/", denominator)), color = "black", size = 6 / .pt, nudge_y = 0.28) +
  facet_wrap(~ trait_class, nrow = 1L) +
  scale_color_manual(values = c(
    replicated_accessible = "#00695C", source_dependent = "#7B1FA2",
    partial = "#518DC9", indeterminate = "#9E9E9E", untestable = "#252525"
  ), drop = FALSE) +
  scale_size_area(max_size = 5, name = "Pairs") +
  labs(
    x = NULL, y = NULL, color = NULL,
    title = "Variant-consistent lineage observability by trait scope"
  ) +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    text = element_text(size = 6, face = "plain"),
    plot.title = element_text(size = 6, face = "plain"),
    axis.text = element_text(size = 6, face = "plain", color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.title = element_text(size = 6, face = "plain"),
    strip.text = element_text(size = 6, face = "plain"),
    strip.background = element_blank(),
    legend.title = element_text(size = 6, face = "plain"),
    legend.text = element_text(size = 6, face = "plain"),
    legend.position = "bottom"
  )
pdf_path <- file.path(OUT, "panels", "fig4d_variant_consistent_lineage_context.pdf")
device <- if (capabilities("cairo")) cairo_pdf else pdf
ggsave(pdf_path, plot, width = 5.4, height = 2.8, device = device)
manifest <- data.table(
  release_id = RELEASE_ID,
  panel = "fig4d_variant_consistent_lineage_context",
  pdf = substring(pdf_path, nchar(CANDIDATE) + 2L),
  source_table = substring(source_path, nchar(CANDIDATE) + 2L),
  pdf_sha256 = digest(pdf_path, algo = "sha256", file = TRUE, serialize = FALSE),
  source_sha256 = digest(source_path, algo = "sha256", file = TRUE, serialize = FALSE)
)
fwrite(manifest, file.path(OUT, "figure_manifest.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
