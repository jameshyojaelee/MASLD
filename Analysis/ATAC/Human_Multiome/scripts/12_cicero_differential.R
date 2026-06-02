#!/usr/bin/env Rscript
# ==========================================================================
# Cicero differential co-accessibility: MASH vs NORMAL
# Runs Cicero separately per condition, then compares
# ==========================================================================

suppressPackageStartupMessages({
  library(cicero)
  library(monocle)
  library(VGAM)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts", "figures", "publication_theme.R"))

CICERO_DIR <- file.path(BASE, "Analysis", "ATAC", "Human_Multiome", "results", "cicero")
FIG_OUT    <- file.path(BASE, "figures", "supplementary", "figS_atac")
peaks      <- readLines(file.path(CICERO_DIR, "peaks.tsv"))

chr_sizes <- data.frame(
  V1 = paste0("chr", c(1:22, "X", "Y")),
  V2 = c(248956422, 242193529, 198295559, 190214555, 181538259,
         170805979, 159345973, 145138636, 138394717, 133797422,
         135086622, 133275309, 114364328, 107043718, 101991189,
         90338345, 83257441, 80373285, 58617616, 64444167,
         46709983, 50818468, 156040895, 57227415)
)

save_panel <- function(p, name, w = fig_half_width, h = 3) {
  path <- file.path(FIG_OUT, name)
  save_fig(p, path, width = w, height = h)
  message("  Saved ", path)
}

cat("==============================================================\n")
cat("Cicero Differential Co-accessibility\n")
cat("==============================================================\n")

# Run Cicero for MASH and NORMAL
cond_conns <- list()
for (cond in c("MASH", "NORMAL")) {
  triplet_file <- file.path(CICERO_DIR, paste0("cicero_input_", cond, ".tsv"))
  message(sprintf("\nProcessing %s...", cond))

  input_c <- fread(triplet_file, header = FALSE,
                   col.names = c("Peak", "Cell", "Count"))
  peaks_c <- unique(input_c$Peak)
  cells_c <- unique(input_c$Cell)
  message(sprintf("  %d cells, %d peaks", length(cells_c), length(peaks_c)))

  mat_c <- sparseMatrix(
    i = match(input_c$Peak, peaks_c),
    j = match(input_c$Cell, cells_c),
    x = rep(1L, nrow(input_c)),
    dims = c(length(peaks_c), length(cells_c)),
    dimnames = list(peaks_c, cells_c)
  )

  fd_c <- new("AnnotatedDataFrame",
              data = data.frame(site_name = peaks_c, row.names = peaks_c))
  pd_c <- new("AnnotatedDataFrame",
              data = data.frame(cells = cells_c, row.names = cells_c))
  cds_c <- newCellDataSet(mat_c, phenoData = pd_c, featureData = fd_c,
                           expressionFamily = VGAM::negbinomial.size())

  cds_c <- detectGenes(cds_c)
  cds_c <- estimateSizeFactors(cds_c)

  set.seed(42)
  cds_c <- reduceDimension(cds_c, max_components = 2, num_dim = 6,
                            reduction_method = "tSNE", norm_method = "none")

  tsne_c <- t(reducedDimA(cds_c))
  colnames(tsne_c) <- c("tSNE_1", "tSNE_2")
  # Ensure rownames match CDS colnames (reduceDimension may reorder)
  rownames(tsne_c) <- colnames(cds_c)

  cicero_c <- make_cicero_cds(cds_c, reduced_coordinates = tsne_c, k = 30)

  message(sprintf("  Running Cicero for %s...", cond))
  conns_c <- run_cicero(cicero_c, chr_sizes, window = 500000,
                        sample_num = 100, silent = TRUE)
  cond_conns[[cond]] <- as.data.table(conns_c)
  message(sprintf("  %s: %d connections", cond, nrow(conns_c)))
}

# Compare MASH vs NORMAL
message("\nComputing differential co-accessibility (MASH - NORMAL)...")
diff <- merge(
  cond_conns[["MASH"]][, .(Peak1, Peak2, coaccess_MASH = coaccess)],
  cond_conns[["NORMAL"]][, .(Peak1, Peak2, coaccess_NORMAL = coaccess)],
  by = c("Peak1", "Peak2"), all = TRUE
)
diff[is.na(coaccess_MASH), coaccess_MASH := 0]
diff[is.na(coaccess_NORMAL), coaccess_NORMAL := 0]
diff[, delta := coaccess_MASH - coaccess_NORMAL]
diff <- diff[order(-abs(delta))]

fwrite(diff, file.path(CICERO_DIR, "differential_coaccessibility.csv"))

n_gained <- sum(diff$delta > 0.1, na.rm = TRUE)
n_lost   <- sum(diff$delta < -0.1, na.rm = TRUE)
message(sprintf("  Total pairs: %s", formatC(nrow(diff), big.mark = ",")))
message(sprintf("  Gained (delta > 0.1): %d", n_gained))
message(sprintf("  Lost (delta < -0.1): %d", n_lost))

# --- Panel 21: Differential co-accessibility ---
message("\nGenerating Panel 21...")
set.seed(42)
p21 <- ggplot(diff[sample(.N, min(.N, 50000))], aes(x = delta)) +
  geom_histogram(bins = 100, fill = masld_colors$down, color = NA, alpha = 0.8) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray40") +
  labs(x = expression(Delta~"co-accessibility (MASH - Normal)"),
       y = "Peak pairs",
       title = "Differential co-accessibility") +
  annotate("text", x = 0.3, y = Inf, vjust = 2,
           label = paste0("Gained: ", n_gained),
           size = 2.2, color = masld_colors$up) +
  annotate("text", x = -0.3, y = Inf, vjust = 2,
           label = paste0("Lost: ", n_lost),
           size = 2.2, color = masld_colors$down) +
  theme_masld()

save_panel(p21, "panel_21_differential_coaccess.pdf", w = fig_half_width, h = 2.5)

cat("\n==============================================================\n")
cat("Differential analysis complete.\n")
cat("==============================================================\n")
