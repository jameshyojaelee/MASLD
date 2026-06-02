#!/usr/bin/env Rscript
#' 307f: Virtual staging v2 — Focused on Hepatocytes + Macrophages only.
#'
#' Shows how bulk fibrosis/NAS transition signatures change along pseudotime
#' for the two cell types where pseudotime is biologically meaningful.
#' Shared y-axis, thicker smoothing, only UP signatures.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime")
FIG_DIR <- file.path(BASE, "figures/supplementary/figS02_progression")

theme_src <- file.path(BASE, "scripts/figures/publication_theme.R")
if (file.exists(theme_src)) source(theme_src)

cat("=== 307f: Virtual Staging v2 ===\n")

# Load and merge
scores <- fread(file.path(PT_DIR, "bulk_signature_scores.csv"))
names(scores)[1] <- "cell"
cpt <- fread(file.path(PT_DIR, "consensus_pseudotime_all.csv"))
names(cpt)[1] <- "cell"

dt <- merge(scores, cpt[, .(cell, consensus_pseudotime)], by = "cell")
dt <- dt[!is.na(consensus_pseudotime)]
cat(nrow(dt), "cells with both scores and pseudotime\n")

# Only Hepatocytes and Macrophages
FOCUS <- c("Hepatocytes", "Macrophages")
dt <- dt[cell_type %in% FOCUS]
cat(nrow(dt), "cells in focus cell types\n")

# Fibrosis UP signatures
fib_cols <- grep("^sig_transition_F.*_UP$", names(dt), value = TRUE)
nas_cols <- grep("^sig_transition_NAS.*_UP$", names(dt), value = TRUE)
cat("Fibrosis UP cols:", fib_cols, "\n")
cat("NAS UP cols:", nas_cols, "\n")

# Bin pseudotime (30 bins for smoother curves)
bin_and_aggregate <- function(sub_dt, sig_cols, n_bins = 30) {
  brks <- unique(quantile(sub_dt$consensus_pseudotime,
                            probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(brks) < 3) return(NULL)
  sub_dt[, pt_bin := as.integer(cut(consensus_pseudotime, breaks = brks,
                                     include.lowest = TRUE, labels = FALSE))]
  # Bin midpoints
  midpoints <- sub_dt[, .(pt_mid = mean(consensus_pseudotime)), by = pt_bin]

  results <- list()
  for (col in sig_cols) {
    binned <- sub_dt[, .(mean_score = mean(get(col), na.rm = TRUE)), by = pt_bin]
    binned <- merge(binned, midpoints, by = "pt_bin")
    binned$signature <- col
    results[[col]] <- binned
  }
  rbindlist(results)
}

# Process each cell type
all_fib <- list()
all_nas <- list()
for (ct in FOCUS) {
  sub <- dt[cell_type == ct]
  cat(ct, ":", nrow(sub), "cells\n")

  fib <- bin_and_aggregate(sub, fib_cols)
  if (!is.null(fib)) { fib$cell_type <- ct; all_fib[[ct]] <- fib }

  nas <- bin_and_aggregate(sub, nas_cols)
  if (!is.null(nas)) { nas$cell_type <- ct; all_nas[[ct]] <- nas }
}

fib_df <- rbindlist(all_fib)
nas_df <- rbindlist(all_nas)

# Clean signature labels
fib_df[, sig_label := gsub("sig_transition_", "", signature)]
fib_df[, sig_label := gsub("_to_", " → ", sig_label)]
fib_df[, sig_label := gsub("_UP", "", sig_label)]
fib_df[, sig_label := factor(sig_label,
                               levels = c("F0 → F1", "F1 → F2", "F2 → F3", "F3 → F4"))]

nas_df[, sig_label := gsub("sig_transition_", "", signature)]
nas_df[, sig_label := gsub("_to_", " → ", sig_label)]
nas_df[, sig_label := gsub("_UP", "", sig_label)]
nas_df[, sig_label := gsub("NAS01", "NAS 0-1", sig_label)]
nas_df[, sig_label := gsub("NAS24", "NAS 2-4", sig_label)]
nas_df[, sig_label := gsub("NAS5", "NAS 5", sig_label)]
nas_df[, sig_label := gsub("NAS68", "NAS 6-8", sig_label)]
nas_df[, sig_label := factor(sig_label,
                               levels = c("NAS 0-1 → NAS 2-4", "NAS 2-4 → NAS 5",
                                           "NAS 5 → NAS 6-8"))]

# Colors: metabolic (blue) → inflammatory (red) progression
fib_colors <- c("F0 → F1" = "#90CAF9", "F1 → F2" = "#42A5F5",
                "F2 → F3" = "#E91E63", "F3 → F4" = "#880E4F")
nas_colors <- c("NAS 0-1 → NAS 2-4" = "#90CAF9",
                "NAS 2-4 → NAS 5" = "#AB47BC",
                "NAS 5 → NAS 6-8" = "#880E4F")

# =========================================================================
# Panel a: Fibrosis signatures along pseudotime
# =========================================================================
pa <- ggplot(fib_df, aes(x = pt_mid, y = mean_score, color = sig_label)) +
  geom_smooth(method = "loess", span = 0.4, se = TRUE, alpha = 0.15, linewidth = 1) +
  scale_color_manual(values = fib_colors, name = "Fibrosis\ntransition") +
  facet_wrap(~cell_type, nrow = 1, scales = "fixed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  theme_masld(base_size = 7) +
  theme(strip.text = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6)) +
  labs(x = NULL, y = "Mean fibrosis\nsignature score",
       title = "Fibrosis program activation along sc pseudotime")

# =========================================================================
# Panel b: NAS signatures along pseudotime
# =========================================================================
pb <- ggplot(nas_df, aes(x = pt_mid, y = mean_score, color = sig_label)) +
  geom_smooth(method = "loess", span = 0.4, se = TRUE, alpha = 0.15, linewidth = 1) +
  scale_color_manual(values = nas_colors, name = "NAS\ntransition") +
  facet_wrap(~cell_type, nrow = 1, scales = "fixed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  theme_masld(base_size = 7) +
  theme(strip.text = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6)) +
  labs(x = "Pseudotime →", y = "Mean NAS\nsignature score",
       title = "NAS program activation along sc pseudotime")

# =========================================================================
# Combine
# =========================================================================
fig <- pa / pb +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

out <- file.path(FIG_DIR, "figS_virtual_staging_v2.pdf")
save_fig(fig, out, width = 5.5, height = 5)
cat("Saved:", out, "\n")

# Print peak ordering for validation
cat("\n--- Peak pseudotime validation ---\n")
for (ct in FOCUS) {
  cat(ct, ":\n")
  peaks <- fib_df[cell_type == ct, .(peak_pt = pt_mid[which.max(mean_score)]),
                   by = sig_label]
  print(peaks[order(peak_pt)])
}

cat("\n=== 307f COMPLETE ===\n")
