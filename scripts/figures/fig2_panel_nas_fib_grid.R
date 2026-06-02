#!/usr/bin/env Rscript
# ============================================================================
# fig2_panel_nas_fib_grid.R
# Fig 2 panel — NAS x Fibrosis staged-sample grid with marginal totals.
# Restored from archive/legacy_pre_redesign_2026-05-17/fig2a.pdf (2026-05-19).
# Panel suffix dropped per convention; Illustrator handles labeling.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(PANEL_DIR, "fig2_panel_nas_fib_grid.pdf")

META_FILE <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
QC_FILE   <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")

BASE_SIZE <- 7
LBL_SIZE  <- 7 / ggplot2::.pt

meta <- fread(META_FILE)
qc   <- fread(QC_FILE)
meta <- merge(meta, qc[, .(sample_id, pass_technical)], by = "sample_id", all.x = TRUE)
meta_qc <- meta[pass_technical == TRUE]

ag <- meta_qc[, .(sample_id, fibrosis_stage, nas_score)]
ag[, fib_int := suppressWarnings(as.integer(fibrosis_stage))]
ag[, nas_int := suppressWarnings(as.integer(nas_score))]
ag_grid <- ag[!is.na(fib_int) & !is.na(nas_int) &
              fib_int %in% 0:4 & nas_int %in% 0:8]
n_both <- nrow(ag_grid)

cell_n    <- ag_grid[, .N, by = .(fib_int, nas_int)]
grid_full <- CJ(fib_int = 0:4, nas_int = 0:8)
grid_full <- merge(grid_full, cell_n, by = c("fib_int", "nas_int"), all.x = TRUE)
grid_full[is.na(N), N := 0]

fib_margin <- grid_full[, .(N = sum(N)), by = fib_int][order(fib_int)]
nas_margin <- grid_full[, .(N = sum(N)), by = nas_int][order(nas_int)]

fib_lvls_ext <- c(paste0("F", 0:4), "Total")
nas_lvls_ext <- c(as.character(0:8), "Total")

main_df   <- grid_full[, .(fib_lab = paste0("F", fib_int),
                            nas_lab = as.character(nas_int),
                            N = N, kind = "cell")]
right_df  <- nas_margin[, .(fib_lab = "Total",
                             nas_lab = as.character(nas_int),
                             N = N, kind = "rowsum")]
top_df    <- fib_margin[, .(fib_lab = paste0("F", fib_int),
                             nas_lab = "Total",
                             N = N, kind = "colsum")]
corner_df <- data.table(fib_lab = "Total", nas_lab = "Total",
                        N = n_both, kind = "grand")

ext <- rbind(main_df, right_df, top_df, corner_df)
ext[, fib_lab := factor(fib_lab, levels = fib_lvls_ext)]
ext[, nas_lab := factor(nas_lab, levels = nas_lvls_ext)]
ext[, fill_n  := ifelse(kind == "cell", N, NA_real_)]

panel <- ggplot(ext, aes(x = fib_lab, y = nas_lab)) +
  geom_tile(aes(fill = fill_n), color = "white", linewidth = 0.4) +
  geom_tile(data = ext[kind != "cell"],
            fill = "#FAFAFA", color = "white", linewidth = 0.4) +
  geom_text(data = ext[kind == "cell"],
            aes(label = ifelse(N == 0, "", as.character(N))),
            size = LBL_SIZE, color = "black") +
  geom_text(data = ext[kind != "cell"],
            aes(label = as.character(N)),
            size = LBL_SIZE, color = "black", fontface = "bold") +
  scale_fill_gradient(low = "#f7fbff", high = "#1565C0",  # deep blue (original)
                      name = "N samples", na.value = "#FAFAFA",
                      guide = guide_colorbar(barwidth = 0.3, barheight = 3)) +
  labs(title = sprintf("Fibrosis x NAS staged samples (n=%s)", comma(n_both)),
       x = "Fibrosis stage", y = "NAS score") +
  theme_masld(base_size = BASE_SIZE) +
  theme(plot.title    = element_text(size = BASE_SIZE),
        axis.title    = element_text(size = BASE_SIZE),
        axis.text     = element_text(size = BASE_SIZE),
        legend.title  = element_text(size = BASE_SIZE),
        legend.text   = element_text(size = BASE_SIZE),
        plot.subtitle = element_blank(),
        plot.margin   = margin(3, 3, 3, 3))

save_fig(panel, OUT_PDF,
         width = fig_half_width * 1.4, height = 3.6)
fwrite(grid_full,  file.path(DATA_DIR, "fig2_panel_nas_fib_grid_data.csv"))
fwrite(fib_margin, file.path(DATA_DIR, "fig2_panel_nas_fib_grid_fib_margin.csv"))
fwrite(nas_margin, file.path(DATA_DIR, "fig2_panel_nas_fib_grid_nas_margin.csv"))
cat(sprintf("[saved] %s\n", OUT_PDF))
