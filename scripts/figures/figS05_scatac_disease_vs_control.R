#!/usr/bin/env Rscript
# figS05_scatac_disease_vs_control.R
#
# Supplementary Fig S4g / Fig S05 diagnostic for human liver snATAC
# differential accessibility. Two panels:
#   (a) one MASLD-vs-control volcano per major cell type from DONOR-LEVEL
#       pseudobulk edgeR-QLF (Squair et al. 2021; effective n = 18 donors,
#       13 MASLD [MASL+MASH] vs 5 NORMAL) — NOT the old per-cell SnapATAC2
#       diff_test (which was pseudoreplicated: cells treated as replicates).
#   (b) one signed bar summary of significant opening/closing peaks across the
#       donor-level MASLD-vs-control contrast and the pseudobulk stage contrasts.
#
# Substrate (rebuilt 2026-07-09): scatac_da_corrected_{hep,stellate,macrophage,
# cholangiocyte}_edger.csv from 07c_corrected_da_celltypes.py + 17c_edger_
# celltypes.R (hep from 17_pseudobulk_de_retest.R). The experimental unit is the
# DONOR throughout, so the significance counts are honest, not cell-count
# artifacts. With only 5 control donors this compartment set remains
# underpowered — report as a donor-level diagnostic, not a headline claim.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrastr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF <- file.path(FIGS05_DIR, "figS05_scatac_disease_vs_control.pdf")
OUT_PDF_FIG4 <- file.path(FIG4_DIR, "panels", "figS4g.pdf")
OUT_CSV <- sub("\\.pdf$", ".summary.csv", OUT_PDF_FIG4)

MAG  <- "#C9265E"
BLUE <- "#1565C0"
GRAY <- "#9E9E9E"

clean_ct <- function(x) {
  map <- c(
    "Hepatocyte" = "Hepatocyte",
    "Stellate_Cell" = "Stellate",
    "Macrophage" = "Macrophage/Kupffer",
    "Kupffer_Cell" = "Macrophage/Kupffer",
    "Cholangiocyte" = "Cholangiocyte",
    "Endothelial" = "Endothelial",
    "LSEC" = "LSEC",
    "B_Cell" = "B cell",
    "NK_T_Cell" = "NK/T cell",
    "Plasma_Cell" = "Plasma cell",
    "Hep" = "Hepatocyte",
    "Fib" = "Stellate",
    "Mac" = "Macrophage/Kupffer",
    "Chol" = "Cholangiocyte",
    "Endo" = "Endothelial"
  )
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

ct_order <- c("Hepatocyte", "Stellate", "Macrophage/Kupffer", "Cholangiocyte")

# ── Panel a: donor-level pseudobulk edgeR-QLF volcanoes per major cell type ──
# One CSV per compartment (feature, logFC, PValue, FDR), effective n = donors.
# Panel a shows the signal-bearing MASH-vs-Control endpoint; the pooled
# MASLD-vs-Control counts are still summarised alongside it in panel b.
SNAP_DIR <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/snapatac2")
masld_files <- c(
  "Hepatocyte"         = "scatac_da_corrected_hep_edger.csv",
  "Stellate"           = "scatac_da_corrected_stellate_edger.csv",
  "Macrophage/Kupffer" = "scatac_da_corrected_macrophage_edger.csv",
  "Cholangiocyte"      = "scatac_da_corrected_cholangiocyte_edger.csv"
)
mash_files <- c(
  "Hepatocyte"         = "scatac_da_mashvsnormal_hep_edger.csv",
  "Stellate"           = "scatac_da_mashvsnormal_stellate_edger.csv",
  "Macrophage/Kupffer" = "scatac_da_mashvsnormal_macrophage_edger.csv",
  "Cholangiocyte"      = "scatac_da_mashvsnormal_cholangiocyte_edger.csv"
)
read_edger <- function(disp, fname) {
  fp <- file.path(SNAP_DIR, fname)
  if (!file.exists(fp)) {
    warning(sprintf("missing donor-level DA file for %s: %s", disp, fp))
    return(NULL)
  }
  d <- fread(fp)
  d[, .(feature, log2FC = as.numeric(logFC),
        pvalue = as.numeric(PValue), padj = as.numeric(FDR),
        cell_type = disp)]
}
read_set <- function(files) {
  x <- rbindlist(Map(read_edger, names(files), unname(files)), use.names = TRUE)
  x[, ct_display := factor(cell_type, levels = ct_order)]
  x <- x[!is.na(ct_display)]
  x[, neglog10_padj := pmin(-log10(pmax(padj, 1e-300)), 50)]
  x[, direction := fcase(
    padj < 0.05 & log2FC > 0, "Opening",
    padj < 0.05 & log2FC < 0, "Closing",
    default = "NS"
  )]
  x[, direction := factor(direction, levels = c("NS", "Opening", "Closing"))]
  setorder(x, direction)
  x[]
}
da_masld <- read_set(masld_files)   # panel b MASLD-vs-Control bars only
da       <- read_set(mash_files)    # panel a volcanoes + panel b MASH bars

volc_counts <- da[, .(
  n_open = sum(direction == "Opening", na.rm = TRUE),
  n_close = sum(direction == "Closing", na.rm = TRUE)
), by = ct_display]
volc_counts[, label := sprintf("%d open / %d close", n_open, n_close)]

xmax <- quantile(abs(da$log2FC[is.finite(da$log2FC)]), 0.995, na.rm = TRUE)
xmax <- max(1, min(3.5, xmax))

p_a <- ggplot(da, aes(x = pmax(pmin(log2FC, xmax), -xmax),
                      y = neglog10_padj, colour = direction)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "gray70", linewidth = 0.22) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "gray70", linewidth = 0.22) +
  rasterise(geom_point(size = 0.26, alpha = 0.45, shape = 16), dpi = 300) +
  geom_text(data = volc_counts, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.1,
            size = 5/.pt, lineheight = 0.86, colour = "grey25") +
  facet_wrap(~ ct_display, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = c(NS = GRAY, Opening = MAG, Closing = BLUE),
                      breaks = c("Opening", "Closing", "NS"),
                      name = NULL) +
  scale_x_continuous(name = expression(log[2]*" FC  (MASH vs Control)"),
                     limits = c(-xmax, xmax), expand = expansion(mult = 0.04)) +
  scale_y_continuous(name = expression(-log[10]*" FDR"),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(tag = "a") +
  theme_masld(base_size = 6) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(1.5, "pt"),
        strip.text = element_text(size = 6, face = "plain",
                                  margin = margin(1, 0, 1, 0)),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6, face = "plain"),
        legend.position = "bottom",
        legend.key.size = unit(0.22, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.5, "pt"),
        legend.text = element_text(size = 6),
        plot.tag = element_text(size = 9, face = "plain"),
        plot.tag.position = c(0.004, 0.996),
        plot.margin = margin(1, 1, 0, 1))

# ── Panel b: signed significant-peak counts for the two DOCUMENTED donor-level
#    contrasts — MASLD vs Control (MASL+MASH, 13 vs 5) and the extreme-severity
#    endpoint MASH vs Control (documented fibrosis 2-4, 9 vs 5). The inferred-
#    F-stage stage_da contrasts (F0/F3/F4) are intentionally NOT shown: F-stage
#    here is scvi_predicted only (0 documented), not a defensible cohort axis. ──
count_dir <- function(x, contrast_label) {
  cc <- x[, .(open  = sum(padj < 0.05 & log2FC > 0, na.rm = TRUE),
              close = sum(padj < 0.05 & log2FC < 0, na.rm = TRUE)), by = ct_display]
  cc[, `:=`(contrast = contrast_label, cell_type = as.character(ct_display))]
  cc[]
}
disease_counts <- count_dir(da_masld, "MASLD vs Control")  # 13 MASLD vs 5 NORMAL
mash_counts    <- count_dir(da,       "MASH vs Control")   # 9 MASH vs 5 NORMAL

bar <- rbind(
  disease_counts[, .(cell_type, contrast, open, close)],
  mash_counts[, .(cell_type, contrast, open, close)],
  fill = TRUE
)
bar[, cell_type := factor(cell_type, levels = ct_order)]
bar <- bar[!is.na(cell_type)]
contrast_order <- c("MASLD vs Control", "MASH vs Control")
bar[, contrast := factor(contrast, levels = contrast_order)]
bar[, total := open + close]

bar_long <- rbindlist(list(
  bar[, .(cell_type, contrast, direction = "Opening", n = open)],
  bar[, .(cell_type, contrast, direction = "Closing", n = -close)]
))
bar_long[, direction := factor(direction, levels = c("Opening", "Closing"))]

p_b <- ggplot(bar_long, aes(x = contrast, y = n, fill = direction)) +
  geom_col(width = 0.78) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
  facet_wrap(~ cell_type, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = c(Opening = MAG, Closing = BLUE), name = NULL) +
  scale_y_continuous(labels = function(x) abs(x), expand = expansion(mult = c(0.14, 0.16))) +
  labs(tag = "b", x = NULL, y = "Significant DA peaks (FDR < 0.05)") +
  theme_masld(base_size = 6) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(1.5, "pt"),
        strip.text = element_text(size = 6, face = "plain",
                                  margin = margin(1, 0, 1, 0)),
        axis.text.y = element_text(size = 5, colour = "black"),
        axis.text.x = element_text(size = 5, angle = 25, hjust = 1, colour = "black"),
        axis.title.y = element_text(size = 6, face = "plain"),
        legend.position = "bottom",
        legend.key.size = unit(0.22, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.5, "pt"),
        legend.text = element_text(size = 6),
        plot.tag = element_text(size = 9, face = "plain"),
        plot.tag.position = c(0.004, 0.996),
        plot.margin = margin(0, 1, 1, 1))

combo <- p_a / p_b + plot_layout(heights = c(1.05, 1.0))
ggsave(OUT_PDF, combo, width = 6.6, height = 4.1, device = cairo_pdf)
ggsave(OUT_PDF_FIG4, combo, width = 6.6, height = 4.1, device = cairo_pdf)
fwrite(bar[, .(cell_type, contrast, open, close, total)], OUT_CSV)
cat(sprintf("[saved] %s\n", OUT_PDF))
cat(sprintf("[saved] %s\n", OUT_PDF_FIG4))
cat(sprintf("[saved] %s\n", OUT_CSV))

message(strrep("=", 78))
message("FIGURE LEGEND (paste into manuscript; stats live here, not on the panel)")
message(strrep("=", 78))
message("Human liver snATAC differential-accessibility diagnostic across major cell types, computed at ",
        "DONOR resolution (donor-level pseudobulk edgeR-QLF; Squair et al. 2021), replacing the earlier ",
        "pseudoreplicated per-cell SnapATAC2 test. Panel a shows MASH-vs-Control volcanoes — the documented ",
        "extreme-severity endpoint (9 MASH, CRN fibrosis stage 2-4, vs 5 NORMAL) — for hepatocyte, stellate, ",
        "macrophage/Kupffer, and cholangiocyte compartments (points rasterised). Panel b summarizes significant ",
        "opening and closing peaks for the two DOCUMENTED contrasts per compartment: the pooled MASLD vs Control ",
        "(MASL+MASH, 13 vs 5) and the MASH-vs-Control endpoint. Dropping the intermediate MASL donors sharpens the ",
        "contrast and recovers the stellate (HSC-activation) signal that the pooled comparison dilutes to zero. A ",
        "literal F0-vs-F4 contrast is not shown because F-stage in this cohort is scVI-inferred only (0 documented) ",
        "and is not a defensible axis. Positive bars are opening peaks, negative bars are closing peaks. Reported as ",
        "a donor-level diagnostic: with only 5 control donors this set is underpowered, so counts index effect size, ",
        "not a headline claim.")
message("\nMASH-vs-Control (panel a) donor-level edgeR counts:")
print(mash_counts[order(ct_display)])
message("\nMASLD-vs-Control donor-level edgeR counts:")
print(disease_counts[order(ct_display)])
message("\nAll-count summary written to: ", OUT_CSV)
message(strrep("=", 78))
