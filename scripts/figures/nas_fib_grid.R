#!/usr/bin/env Rscript
# ============================================================================
# nas_fib_grid.R
# Fig 3 panel (figs3b) — Fibrosis x NAS staged-sample grid with marginal totals.
# Square cells (coord_fixed), gapped marginal band, contrast-aware labels.
# Restored from archive/legacy_pre_redesign_2026-05-17/fig2a.pdf (2026-05-19).
# Output: FIG2_DIR/panels/figs3b_nas_fib_grid.pdf
#   (FIG2_DIR resolves to figures/main/fig3_RNAseq — known back-compat misnomer)
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

PANEL_DIR <- file.path(FIG2_DIR, "panels")   # main fig3_RNAseq
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)
OUT_PDF <- file.path(PANEL_DIR, "figs3b_nas_fib_grid.pdf")

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

# ── Numeric tile positions with a small gap before the marginal "Total" band ──
GAP <- 0.45
ext[, x := as.numeric(fib_lab)]            # F0..F4 = 1..5 ; Total = 6
ext[, y := as.numeric(nas_lab)]            # 0..8    = 1..9 ; Total = 10
ext[fib_lab == "Total", x := 6 + GAP]      # push Total column right of the grid
ext[nas_lab == "Total", y := 10 + GAP]     # push Total row above the grid
x_breaks <- c(1:5, 6 + GAP)
y_breaks <- c(1:9, 10 + GAP)

# Contrast-aware cell labels: white on the dark-blue (high-count) cells.
ext[, txt_col := "grey15"]
ext[kind == "cell" & N >= 28, txt_col := "white"]

# ── Plot ─────────────────────────────────────────────────────────────────────
panel <- ggplot(ext, aes(x = x, y = y)) +
  # marginal Total band — neutral fill, drawn first
  geom_tile(data = ext[kind != "cell"],
            width = 1, height = 1,
            fill = "#ECEFF1", color = "white", linewidth = 0.6) +
  # main grid — blue gradient by sample count
  geom_tile(data = ext[kind == "cell"],
            aes(fill = fill_n),
            width = 1, height = 1, color = "white", linewidth = 0.6) +
  # cell counts (blank for zero)
  geom_text(data = ext[kind == "cell" & N > 0],
            aes(label = N, color = txt_col), size = LBL_SIZE) +
  # marginal totals — bold
  geom_text(data = ext[kind != "cell"],
            aes(label = N), color = "grey15",
            fontface = "bold", size = LBL_SIZE) +
  scale_color_identity() +
  scale_fill_gradient(low = "#EAF2FB", high = "#1565C0",
                      name = "Samples (n)", breaks = pretty_breaks(4),
                      guide = guide_colorbar(barwidth = 0.4, barheight = 4.2,
                                             ticks.colour = "white")) +
  scale_x_continuous(breaks = x_breaks, labels = fib_lvls_ext,
                     expand = expansion(add = 0.06)) +
  scale_y_continuous(breaks = y_breaks, labels = nas_lvls_ext,
                     expand = expansion(add = 0.06)) +
  coord_fixed() +
  labs(title    = "Sample staging composition",
       subtitle = sprintf("Fibrosis stage × NAS score · n = %s paired-staged samples",
                          comma(n_both)),
       x = "Fibrosis stage", y = "NAS score") +
  theme_masld(base_size = BASE_SIZE) +
  theme(panel.grid    = element_blank(),
        panel.border  = element_blank(),
        axis.ticks    = element_blank(),
        plot.title    = element_text(size = BASE_SIZE + 1, face = "bold"),
        plot.subtitle = element_text(size = BASE_SIZE - 0.5, colour = "grey35",
                                     margin = margin(b = 6)),
        axis.title    = element_text(size = BASE_SIZE),
        axis.text     = element_text(size = BASE_SIZE),
        legend.title  = element_text(size = BASE_SIZE - 0.5),
        legend.text   = element_text(size = BASE_SIZE - 0.5),
        plot.margin   = margin(4, 4, 4, 4))

save_fig(panel, OUT_PDF, width = 4.7, height = 6.0)
fwrite(grid_full,  file.path(DATA_DIR, "nas_fib_grid_data.csv"))
fwrite(fib_margin, file.path(DATA_DIR, "nas_fib_grid_fib_margin.csv"))
fwrite(nas_margin, file.path(DATA_DIR, "nas_fib_grid_nas_margin.csv"))
cat(sprintf("[saved] %s\n", OUT_PDF))
