##############################################################################
# Fig 4 (Team A): Mitochondrial spatial condensation in steatosis
#
# Message: mtDNA-OXPHOS transcripts GAIN spatial autocorrelation in steatosis
#   (MT-CO3 Moran's I 0.57 -> 0.81), an unbiased spatial-transcriptomics readout
#   of mitochondrial spatial condensation (peridroplet mitochondria; Benador 2018).
#
# Panel: paired-slope plot of squidpy per-donor Moran's I, Healthy -> Steatotic,
#   for the 12 mtDNA protein-coding transcripts present in differential_svgs.csv.
#
# DATA PROVENANCE (hard-gate compliant — ledger row mito_clustering):
#   Analysis/Spatial/results/svg/differential_svgs.csv
#   columns: morans_I_healthy, morans_I_masld, delta_I
#   metric : squidpy per-donor Moran's I (NOT Hotspot C/Z autocorrelation)
#
# CAVEAT (baked into legend + stdout): squidpy per-donor Moran within GSE192741
#   only (2 Healthy JBO vs 3 Steatotic JBO, n=5); the biology is CONFIRMATORY
#   (peridroplet mitochondria, Benador 2018), the ST metric is the novel part.
#
# HONEST READ: the file shows 11 of 12 transcripts rise; MT-ATP8 is the lone
#   exception (delta_I = -0.081). The "all 12 rise" framing in the prose/ledger
#   is softened here to "11 of 12 (MT-ATP8 excepted)" — we plot what the file says.
#
# Output: figures/main/fig4_validation/mito_spatial_clustering.pdf
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ---------------------------------------------------------------------------
# Load — canonical squidpy SVG table (ledger: mito_clustering)
# ---------------------------------------------------------------------------
svg_path <- file.path(BASE, "Analysis/Spatial/results/svg/differential_svgs.csv")
stopifnot(file.exists(svg_path))
svg <- fread(svg_path)
# first column is the unnamed gene index
setnames(svg, 1, "gene")

# mtDNA protein-coding transcripts only: "MT-" prefix (excludes nuclear MT1*/MTRNR*)
mt <- svg[grepl("^MT-", gene)]
mt <- mt[order(-delta_I)]

stopifnot(nrow(mt) == 12L)

# Honest tally
n_up   <- sum(mt$delta_I > 0)
n_down <- sum(mt$delta_I <= 0)
co3    <- mt[gene == "MT-CO3"]
message(sprintf("[mito_clustering] %d mtDNA transcripts read; %d rise, %d fall (MT-ATP8 delta_I=%.3f).",
                nrow(mt), n_up, n_down, mt[gene == "MT-ATP8", delta_I]))
message(sprintf("[mito_clustering] Headline MT-CO3: Moran's I %.3f (Healthy) -> %.3f (Steatotic), delta_I = +%.3f.",
                co3$morans_I_healthy, co3$morans_I_masld, co3$delta_I))

# ---------------------------------------------------------------------------
# Reshape to paired long form
# ---------------------------------------------------------------------------
long <- rbind(
  data.table(gene = mt$gene, cond = "Healthy",   I = mt$morans_I_healthy, delta = mt$delta_I),
  data.table(gene = mt$gene, cond = "Steatotic", I = mt$morans_I_masld,   delta = mt$delta_I)
)
long[, cond := factor(cond, levels = c("Healthy", "Steatotic"))]
# rises = disease-up magenta; the single faller (MT-ATP8) = neutral gray
long[, rises := delta > 0]

col_rise <- masld_colors$up      # #C9265E
col_fall <- masld_colors$ns      # #9E9E9E  (control/neutral invariant)
col_hi   <- masld_colors$mash    # same magenta family, headline gene

# right-edge labels at the Steatotic x-position
lab <- long[cond == "Steatotic"]

# ---------------------------------------------------------------------------
# Plot — paired slopes
# ---------------------------------------------------------------------------
p <- ggplot(long, aes(x = cond, y = I, group = gene)) +
  geom_line(aes(color = rises, linewidth = gene == "MT-CO3"), alpha = 0.85) +
  geom_point(aes(color = rises), size = 1.1) +
  ggrepel::geom_text_repel(
    data = lab,
    aes(label = gene, color = rises),
    size = PUB_GEOM_TEXT, hjust = 0, direction = "y", nudge_x = 0.08,
    segment.size = 0.2, segment.color = "gray70", segment.alpha = 0.6,
    box.padding = 0.08, min.segment.length = 0, max.overlaps = Inf,
    xlim = c(2.05, NA), show.legend = FALSE
  ) +
  scale_color_manual(values = c(`TRUE` = col_rise, `FALSE` = col_fall), guide = "none") +
  scale_linewidth_manual(values = c(`TRUE` = 1.0, `FALSE` = 0.45), guide = "none") +
  scale_x_discrete(expand = expansion(mult = c(0.10, 0.42))) +
  scale_y_continuous(limits = c(0, 0.92), breaks = seq(0, 0.9, 0.2),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Moran's I (spatial autocorrelation)") +
  theme_masld() + theme_pub() +
  theme(
    axis.text.x = element_text(size = PUB_AXIS_TITLE),
    plot.margin = margin(4, 4, 4, 4)
  )

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_pdf <- file.path(FIG4_DIR, "_supp", "mito_spatial_clustering.pdf")
dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
# useDingbats=FALSE set globally via pdf.options() in publication_theme.R
ggsave(out_pdf, p, width = 55 / 25.4, height = 62 / 25.4, device = cairo_pdf)
message("[mito_clustering] wrote ", out_pdf)

# ---------------------------------------------------------------------------
# Legend text (stats live here, not on the panel)
# ---------------------------------------------------------------------------
message("\n==== LEGEND (mito_spatial_clustering) ====")
message(sprintf(
  paste0("Mitochondrial spatial condensation in steatosis. Paired squidpy per-donor ",
         "Moran's I for the 12 mtDNA protein-coding transcripts in the Visium SVG ",
         "table, Healthy -> Steatotic. %d of 12 transcripts GAIN spatial ",
         "autocorrelation (magenta; MT-ATP8 the lone exception, gray, delta_I = %.3f); ",
         "MT-CO3 rises from %.2f to %.2f (delta_I = +%.2f, bold). The coordinated ",
         "rise is an unbiased spatial-transcriptomics readout of mitochondrial ",
         "spatial condensation (peridroplet mitochondria; Benador et al. 2018). ",
         "CAVEAT: squidpy per-donor Moran (not Hotspot C/Z autocorrelation), within ",
         "GSE192741 only (2 Healthy vs 3 Steatotic donors, n=5); the biology is ",
         "confirmatory, the spatial-transcriptomics metric is the novel contribution."),
  n_up, mt[gene == "MT-ATP8", delta_I],
  co3$morans_I_healthy, co3$morans_I_masld, co3$delta_I))
message("==========================================\n")
