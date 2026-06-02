#!/usr/bin/env Rscript
#' 307e: Virtual Staging — bulk fibrosis/NAS transition signatures along
#'       single-cell pseudotime.
#'
#' Key biological insight: as cells progress along pseudotime, which disease
#' programs activate and in what order?  If F0->F1 UP peaks early and F3->F4 UP
#' peaks late, pseudotime recapitulates fibrosis progression at single-cell
#' resolution.
#'
#' Inputs
#'   - bulk_signature_scores.csv   (305 output: per-cell transition scores)
#'   - consensus_pseudotime_all.csv (consensus pseudotime per cell)
#'
#' Output
#'   - figures/supplementary/figS02_progression/figS_virtual_staging.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PT_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime")
FIG_DIR <- file.path(BASE, "figures/supplementary/figS02_progression")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

theme_src <- file.path(BASE, "scripts/figures/publication_theme.R")
if (file.exists(theme_src)) source(theme_src)

# =========================================================================
# Constants
# =========================================================================
CELL_TYPES <- c("Hepatocytes", "Macrophages", "Fibroblasts",
                "Endothelial_cells", "Cholangiocytes")
CT_LABELS  <- c(Hepatocytes = "Hepatocytes", Macrophages = "Macrophages",
                Fibroblasts = "Fibroblasts", Endothelial_cells = "Endothelial",
                Cholangiocytes = "Cholangiocytes")

N_BINS <- 20

# Fibrosis UP signatures (light blue -> dark magenta: metabolic-to-inflammatory)
FIB_SIGS <- c(
  "sig_transition_F0_to_F1_UP",
  "sig_transition_F1_to_F2_UP",
  "sig_transition_F2_to_F3_UP",
  "sig_transition_F3_to_F4_UP"
)
FIB_LABELS <- c(
  sig_transition_F0_to_F1_UP = "F0 \u2192 F1",
  sig_transition_F1_to_F2_UP = "F1 \u2192 F2",
  sig_transition_F2_to_F3_UP = "F2 \u2192 F3",
  sig_transition_F3_to_F4_UP = "F3 \u2192 F4"
)
FIB_COLORS <- c(
  "F0 \u2192 F1" = "#90CAF9",
  "F1 \u2192 F2" = "#42A5F5",
  "F2 \u2192 F3" = "#E91E63",
  "F3 \u2192 F4" = "#880E4F"
)

# NAS UP signatures (analogous gradient)
NAS_SIGS <- c(
  "sig_transition_NAS01_to_NAS24_UP",
  "sig_transition_NAS24_to_NAS5_UP",
  "sig_transition_NAS5_to_NAS68_UP"
)
NAS_LABELS <- c(
  sig_transition_NAS01_to_NAS24_UP = "NAS 0\u20131 \u2192 2\u20134",
  sig_transition_NAS24_to_NAS5_UP  = "NAS 2\u20134 \u2192 5",
  sig_transition_NAS5_to_NAS68_UP  = "NAS 5 \u2192 6\u20138"
)
NAS_COLORS <- c(
  "NAS 0\u20131 \u2192 2\u20134" = "#90CAF9",
  "NAS 2\u20134 \u2192 5"        = "#AB47BC",
  "NAS 5 \u2192 6\u20138"        = "#880E4F"
)

cat("=== 307e: Virtual Staging ===\n")

# =========================================================================
# Load & merge
# =========================================================================
cat("Loading data...\n")
scores <- fread(file.path(PT_DIR, "bulk_signature_scores.csv"))
names(scores)[1] <- "cell"

pt <- fread(file.path(PT_DIR, "consensus_pseudotime_all.csv"))
names(pt)[1] <- "cell"

# Merge on cell barcode — keep only cells present in both
dat <- merge(scores[, c("cell", "cell_type", FIB_SIGS, NAS_SIGS), with = FALSE],
             pt[, .(cell, consensus_pseudotime)],
             by = "cell")

# Drop cells without pseudotime
dat <- dat[!is.na(consensus_pseudotime)]

cat(sprintf("  %s cells with both scores and pseudotime\n", format(nrow(dat), big.mark = ",")))

# =========================================================================
# Helper: bin, average, smooth, and return long-format table
# =========================================================================
bin_and_smooth <- function(dt, sig_cols, sig_labels, ct_col = "cell_type",
                           pt_col = "consensus_pseudotime", n_bins = N_BINS,
                           span = 0.3) {
  results <- list()
  for (ct in CELL_TYPES) {
    sub <- dt[get(ct_col) == ct]
    if (nrow(sub) < n_bins * 5) next  # skip if too few cells

    # Equal-frequency bins (use unique breaks to handle ties)
    brks <- unique(quantile(sub[[pt_col]], probs = seq(0, 1, length.out = n_bins + 1),
                             na.rm = TRUE))
    if (length(brks) < 3) next  # skip if pseudotime is essentially constant
    sub[, pt_bin := as.integer(cut(get(pt_col), breaks = brks,
                                   include.lowest = TRUE, labels = FALSE))]

    # Bin midpoints (mean pseudotime per bin)
    bin_mid <- sub[, .(pt_mid = mean(get(pt_col))), by = pt_bin]

    for (sig in sig_cols) {
      bin_mean <- sub[, .(mean_score = mean(get(sig), na.rm = TRUE)), by = pt_bin]
      bin_mean <- merge(bin_mean, bin_mid, by = "pt_bin")

      # Loess smooth
      if (nrow(bin_mean) >= 4) {
        lo <- tryCatch(
          loess(mean_score ~ pt_mid, data = bin_mean, span = span),
          error = function(e) NULL
        )
        if (!is.null(lo)) {
          bin_mean[, smoothed := predict(lo, newdata = pt_mid)]
        } else {
          bin_mean[, smoothed := mean_score]
        }
      } else {
        bin_mean[, smoothed := mean_score]
      }

      bin_mean[, `:=`(cell_type = ct,
                      cell_type_label = CT_LABELS[ct],
                      signature = sig_labels[sig])]
      results[[paste0(ct, "_", sig)]] <- bin_mean
    }
  }
  rbindlist(results)
}

# =========================================================================
# Compute binned, smoothed scores
# =========================================================================
cat("Binning and smoothing fibrosis signatures...\n")
fib_df <- bin_and_smooth(dat, FIB_SIGS, FIB_LABELS)
fib_df[, signature := factor(signature, levels = FIB_LABELS)]

cat("Binning and smoothing NAS signatures...\n")
nas_df <- bin_and_smooth(dat, NAS_SIGS, NAS_LABELS)
nas_df[, signature := factor(signature, levels = NAS_LABELS)]

# Order cell type facets
ct_order <- CT_LABELS[CELL_TYPES]
fib_df[, cell_type_label := factor(cell_type_label, levels = ct_order)]
nas_df[, cell_type_label := factor(cell_type_label, levels = ct_order)]

# =========================================================================
# Panel A: Fibrosis transition signatures along pseudotime
# =========================================================================
cat("Building fibrosis panel...\n")
p_fib <- ggplot(fib_df, aes(x = pt_mid, y = smoothed, color = signature)) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(y = mean_score), size = 0.4, alpha = 0.35) +
  facet_wrap(~ cell_type_label, nrow = 1, scales = "free_y") +
  scale_color_manual(values = FIB_COLORS, name = "Fibrosis transition") +
  labs(x = "Pseudotime", y = "Mean UP-signature score") +
  theme_masld(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"))

# =========================================================================
# Panel B: NAS transition signatures along pseudotime
# =========================================================================
cat("Building NAS panel...\n")
p_nas <- ggplot(nas_df, aes(x = pt_mid, y = smoothed, color = signature)) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(y = mean_score), size = 0.4, alpha = 0.35) +
  facet_wrap(~ cell_type_label, nrow = 1, scales = "free_y") +
  scale_color_manual(values = NAS_COLORS, name = "NAS transition") +
  labs(x = "Pseudotime", y = "Mean UP-signature score") +
  theme_masld(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"))

# =========================================================================
# Combine & save
# =========================================================================
cat("Combining panels...\n")
combined <- p_fib / p_nas +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

out_path <- file.path(FIG_DIR, "figS_virtual_staging.pdf")
save_fig(combined, out_path, width = 7.09, height = 5)
cat(sprintf("Saved: %s\n", out_path))

# =========================================================================
# Summary statistics
# =========================================================================
cat("\n--- Fibrosis signature peak pseudotime by cell type ---\n")
fib_peaks <- fib_df[, .(peak_pt = pt_mid[which.max(smoothed)],
                         peak_score = max(smoothed, na.rm = TRUE)),
                    by = .(cell_type_label, signature)]
print(fib_peaks[order(cell_type_label, signature)])

cat("\n--- NAS signature peak pseudotime by cell type ---\n")
nas_peaks <- nas_df[, .(peak_pt = pt_mid[which.max(smoothed)],
                         peak_score = max(smoothed, na.rm = TRUE)),
                    by = .(cell_type_label, signature)]
print(nas_peaks[order(cell_type_label, signature)])

# Ordering test: for each cell type, is the peak of F0->F1 earlier than F3->F4?
cat("\n--- Ordering validation (F0->F1 peak < F3->F4 peak?) ---\n")
for (ct in ct_order) {
  early <- fib_peaks[cell_type_label == ct & signature == "F0 \u2192 F1", peak_pt]
  late  <- fib_peaks[cell_type_label == ct & signature == "F3 \u2192 F4", peak_pt]
  if (length(early) > 0 && length(late) > 0) {
    cat(sprintf("  %s: F0->F1 peak=%.3f, F3->F4 peak=%.3f, ordered=%s\n",
                ct, early, late, ifelse(early < late, "YES", "no")))
  }
}

cat("\n=== 307e complete ===\n")
