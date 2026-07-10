# Fig 3 panel — chromVAR hepatocyte TF activity across MASLD progression
# Same 18 curated TFs as the TF × cell-type panel, but restricted to
# hepatocyte cells and stratified by disease condition (NORMAL / MASL / MASH).
# Donor exclusions match the SCENIC+ donor-regulon panel.
#
# SIGNIFICANCE CAVEAT (honest framing): the point-size -log10(padj) is a
# Wilcoxon-vs-F0 BH-corrected WITHIN this curated ~18-TF hypothesis set — it is
# EXPLORATORY, hypothesis-driven significance, NOT genome-wide. Across the full
# donor-level chromVAR test space (~7,133 TF x cell-type tests) NOTHING survives
# BH correction (n=18 multiome is underpowered; the inflated per-cell "4,832"
# was pseudoreplication — see memory/project-megareview-2026-06-14). Treat this
# panel as a descriptive deviation map on pre-selected TFs, not a discovery claim.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
PER_DONOR_TSV <- file.path(ATAC_DIR, "results/chromvar_v2/chromvar_per_donor.tsv.gz")
META_FILE     <- file.path(ATAC_DIR, "metadata/donor_metadata_curated.tsv")
OUT_PDF       <- file.path(FIGS05_DIR, "figS05_scatac_chromvar_hepatocyte_progression.pdf")
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Curated TF list (display name -> chromVAR motif name) — matches the
# fig3_chromvar_tfct_dotplot panel.
# -----------------------------------------------------------------------------
TF_GROUPS <- list(
  `Hepatocyte-related` =
    c(HNF4A  = "HNF4A",  HNF1A = "HNF1A",
      FOXA1  = "FOXA1",  FOXA2 = "FOXA2"),
  `NR drug targets` =
    c(THRB   = "THRB",   NR1H4 = "Nr1H4",
      PPARA  = "Ppara",  PPARG = "PPARG",
      RORA   = "RORA",   RXRA  = "Rxra"),
  `Metabolic` =
    c(MLXIPL = "MLXIPL", `NR1H3 (LXRα)` = "Nr1h3", KLF15 = "KLF15"),
  `Stress / inflammation` =
    c(XBP1   = "XBP1",   ETS1  = "ETS1"),
  `Fibrogenic / EMT` =
    c(SNAI1  = "SNAI1",  ZEB1  = "ZEB1",  TCF4  = "TCF4")
)

BOLD_TFS <- c("THRB", "NR1H4", "RORA", "HNF4A", "PPARA", "PPARG")

# Donor exclusions (matching the SCENIC+ donor-regulon panel)
DROP_LOW_HEP  <- c("D05", "D07", "D18")   # < 100 hepatocytes in SCENIC+ matrix
EXCLUDE_OUTL  <- c("D15")                 # regulon-preserved cirrhotic

# Disease-condition binning (NORMAL / MASL / MASH)
STAGE_LEVELS  <- c("NORMAL", "MASL", "MASH")

# -----------------------------------------------------------------------------
# Load + filter
# -----------------------------------------------------------------------------
tf_lookup <- rbindlist(lapply(names(TF_GROUPS), function(g) {
  v <- TF_GROUPS[[g]]
  data.table(display = names(v), motif = unname(v), group = g)
}))

pd <- fread(cmd = sprintf("zcat %s", PER_DONOR_TSV))
hep <- pd[cell_type == "Hepatocytes"]
# The per-donor TSV already uses D-codes (e.g., "D01") in donor_id_atac
setnames(hep, "donor_id_atac", "donor_id")

# Bring in F-stage + condition via D-code merge
meta <- fread(META_FILE)
# Backwards-compat for protocol contamination remediation: drop donors flagged
# by `exclude_stage_analysis` (added to donor_metadata_extended.tsv; absent
# from the ATAC-only donor_metadata_curated.tsv used here, hence the fallback).
if (!"exclude_stage_analysis" %in% names(meta)) meta[, exclude_stage_analysis := FALSE]
meta <- meta[exclude_stage_analysis != TRUE]
# Prefer `F_stage_augmented_clean` (post-contamination fix) when available;
# fall back to legacy `F_stage_augmented` for the ATAC curated metadata.
if (!"F_stage_augmented_clean" %in% names(meta)) {
  meta[, F_stage_augmented_clean := F_stage_augmented]
}
hep <- merge(hep, meta[, .(donor_id, condition,
                            F_stage_augmented = F_stage_augmented_clean)],
             by = "donor_id", all.x = FALSE)

cat("Hepatocyte rows:", nrow(hep), "; donors:", uniqueN(hep$donor_id), "\n")

# Apply donor exclusions
drop_set <- c(DROP_LOW_HEP, EXCLUDE_OUTL)
hep <- hep[!donor_id %in% drop_set]
cat("After exclusions (", paste(drop_set, collapse = ","), "):",
    uniqueN(hep$donor_id), "donors\n")

# Bin by MASLD condition
hep[, F_bin := factor(condition, levels = STAGE_LEVELS)]
cat("\nDonors per condition bin:\n")
print(unique(hep[, .(donor_id, condition, F_bin)])[, .N, by = F_bin])

# Restrict to the 18 curated TFs (motif-name match, then collapse motif
# variants — e.g., HNF4A vs HNF4A.1 — by TF_symbol mean)
hep <- hep[TF %in% tf_lookup$motif | TF_symbol %in% tf_lookup$motif]
hep <- hep[, .(mean_deviation = mean(mean_deviation, na.rm = TRUE),
               n_cells        = max(n_cells)),
           by = .(donor_id, TF_symbol, condition, F_stage_augmented, F_bin)]
hep <- merge(hep, tf_lookup, by.x = "TF_symbol", by.y = "motif", all.x = TRUE)
hep <- hep[!is.na(display)]
cat("\nTFs retained:", uniqueN(hep$display), "\n")
print(unique(hep$display))

# -----------------------------------------------------------------------------
# Per-bin summary: mean deviation, Wilcoxon vs F0, n donors
# -----------------------------------------------------------------------------
bin_stats <- hep[, .(
  mean_dev = mean(mean_deviation, na.rm = TRUE),
  n_donors = uniqueN(donor_id),
  values   = list(mean_deviation)
), by = .(display, group, F_bin)]

# Wilcoxon vs F0 (per TF)
ref_vals <- bin_stats[F_bin == STAGE_LEVELS[1], .(display, ref_vals = values)]
bin_stats <- merge(bin_stats, ref_vals, by = "display", all.x = TRUE)
bin_stats[, pval := mapply(function(v, ref) {
  if (length(v) < 2 || length(ref) < 2) return(NA_real_)
  tryCatch(wilcox.test(unlist(v), unlist(ref))$p.value,
           error = function(e) NA_real_)
}, values, ref_vals)]
bin_stats[F_bin == STAGE_LEVELS[1], pval := NA_real_]   # no self-comparison

# BH adjust across the non-reference rows
bin_stats[!is.na(pval), padj := p.adjust(pval, method = "BH")]
bin_stats[, neglog10_padj := pmin(-log10(pmax(padj, 1e-10)), 4)]
bin_stats[F_bin == STAGE_LEVELS[1], neglog10_padj := NA_real_]
bin_stats[, mean_dev_clip := pmax(pmin(mean_dev, 0.6), -0.6)]

# Axis ordering
y_order <- unlist(lapply(TF_GROUPS, names), use.names = FALSE)
y_order <- intersect(y_order, bin_stats$display)
bin_stats[, display := factor(display, levels = rev(y_order))]
bin_stats[, F_bin   := factor(F_bin, levels = STAGE_LEVELS)]

# Bold-face y labels for user-specified TFs
y_levels <- levels(bin_stats$display)
y_stem   <- sub(" \\(.*", "", y_levels)
y_face   <- ifelse(y_stem %in% BOLD_TFS, "italic", "plain")

# Group separators
group_of <- tf_lookup$group[match(y_levels, tf_lookup$display)]
breaks_between <- which(diff(as.integer(factor(group_of, levels = unique(group_of)))) != 0)
hline_y <- length(y_levels) - breaks_between + 0.5

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
p <- ggplot(bin_stats,
            aes(x = F_bin, y = display,
                colour = mean_dev_clip, size = neglog10_padj)) +
  geom_hline(yintercept = hline_y,
             colour = "gray85", linewidth = 0.3) +
  geom_point(alpha = 0.95) +
  # Reference column (NORMAL): plain dots, no size encoding
  geom_point(data = bin_stats[F_bin == STAGE_LEVELS[1]],
             aes(x = F_bin, y = display, colour = mean_dev_clip),
             size = 2.0, inherit.aes = FALSE) +
  scale_colour_gradient2(
    low      = "#1565C0",
    mid      = "#9E9E9E",
    high     = "#C9265E",
    midpoint = 0,
    limits   = c(-0.6, 0.6),
    oob      = scales::squish,
    name     = "Motif deviation z"
  ) +
  scale_size_continuous(
    range  = c(0.6, 3.6),
    limits = c(0, 4),
    breaks = c(1, 2, 3),
    labels = c("1", "2", expression("">=3)),
    name   = expression(-log[10]*" padj")
  ) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(
    axis.text.x      = element_text(size = PUB_AXIS_TEXT, colour = "black"),
    axis.text.y      = element_text(size = PUB_AXIS_TEXT, face = y_face,
                                    colour = "black"),
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    legend.margin    = margin(l = 4, r = 0),
    legend.title     = element_text(size = PUB_LEGEND_TIT, face = "plain"),
    legend.text      = element_text(size = PUB_LEGEND),
    plot.margin      = margin(4, 6, 2, 4)
  ) +
  guides(
    colour = guide_colourbar(title.position = "top",
                             barwidth  = unit(0.22, "cm"),
                             barheight = unit(2.0, "cm"),
                             order = 1),
    size   = guide_legend(title.position = "top",
                          override.aes = list(colour = "#444444"),
                          order = 2)
  )

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
ggsave(OUT_PDF, p,
       width = 3.6, height = 3.6,
       device = cairo_pdf)

# -----------------------------------------------------------------------------
# Source CSV
# -----------------------------------------------------------------------------
out_csv <- sub("\\.pdf$", ".csv", OUT_PDF)
fwrite(bin_stats[, .(display, group, F_bin, mean_dev, n_donors, pval, padj)],
       out_csv)
message(sprintf("[fig3_chromvar_hepatocyte_progression] wrote %s", OUT_PDF))
message(sprintf("  TFs in plot      : %d", uniqueN(bin_stats$display)))
message(sprintf("  F-stage bins     : %s", paste(STAGE_LEVELS, collapse = ", ")))
message(sprintf("  Donors retained  : %d", uniqueN(hep$donor_id)))
message(sprintf("  Source CSV       : %s", out_csv))
