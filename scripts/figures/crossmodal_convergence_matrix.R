#!/usr/bin/env Rscript
# ============================================================================
# crossmodal_convergence_matrix.R
# Main Fig 3 (RNA-seq) — cross-modal convergence dot-heatmap
#
# A rare convergent core: genes that are simultaneously a robust bulk DEG
# (bulk_lfsr<0.05 & |bulk_shrunk_logFC|>0.5) AND a strong hepatic colocalization
# hit (coloc_best_pp4>0.5). For each convergent gene we show which orthogonal
# modalities are "lit": human bulk, mouse bulk, sc-hepatocyte, Hotspot module
# membership, spatial GeoMx direction, and colocalization presence.
#
# INTEGRITY: COLOC PP.H4 is read ONLY from gene_level_coloc.csv (coloc_best_pp4,
# valid <=1). The atlas coloc_susie_best_pp4 column is CORRUPTED (values >1) and
# is NEVER used. A regression assert enforces max(PP.H4) <= 1.
#
# Effect modalities (bulk / mouse / sc-hep): color = signed direction
#   (magenta up / blue down), dot size = |effect|.
# Direction/presence modalities (spatial GeoMx direction, Hotspot, COLOC):
#   rendered as direction or presence only (GeoMx padj all >0.85 -> not a
#   significance claim; Hotspot/COLOC are membership/presence).
#
# Output: figures/main/fig3_RNAseq/panels/fig3i_crossmodal_convergence_matrix.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
OUT_PDF   <- file.path(PANEL_DIR, "fig3i_crossmodal_convergence_matrix.pdf")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
coloc <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))

# INTEGRITY regression check
stopifnot(max(coloc$coloc_best_pp4, na.rm = TRUE) <= 1)

acols <- c("human_symbol", "bulk_shrunk_logFC", "bulk_lfsr", "mouse_meta_logFC",
           "n_diets_sig", "sc_hepatocyte_logFC", "hotspot_n_modules")
m <- merge(atlas[, ..acols], coloc[, .(human_symbol = gene, coloc_best_pp4)],
           by = "human_symbol")

# ---------------------------------------------------------------------------
# Convergent core: robust bulk DEG AND strong hepatic COLOC
# ---------------------------------------------------------------------------
core <- m[bulk_lfsr < 0.05 & abs(bulk_shrunk_logFC) > 0.5 & coloc_best_pp4 > 0.5]

# Per-modality "lit" flags
core[, lit_bulk    := TRUE]                                   # core defn
core[, lit_mouse   := !is.na(mouse_meta_logFC) & (n_diets_sig > 0 |
                        abs(mouse_meta_logFC) > 0.5)]
core[, lit_schep   := !is.na(sc_hepatocyte_logFC) & abs(sc_hepatocyte_logFC) > 0.25]
core[, lit_hotspot := !is.na(hotspot_n_modules) & hotspot_n_modules > 0]
core[, lit_coloc   := TRUE]                                   # core defn

core[, n_lit := lit_bulk + lit_mouse + lit_schep + lit_hotspot +
                lit_coloc]
setorder(core, -n_lit, -coloc_best_pp4)

# ---------------------------------------------------------------------------
# Long format for dot-heatmap (only lit cells get a dot)
# ---------------------------------------------------------------------------
mod_levels <- c("Human bulk", "Mouse bulk", "sc-Hepatocyte",
                "Hotspot", "Coloc PP.H4")

mk <- function(g, mod, lit, eff = NA_real_, kind) {
  if (!lit) return(NULL)
  dir <- if (kind == "effect" || kind == "direction") {
           if (eff > 0) "Up" else "Down"
         } else "Present"
  size <- if (kind == "effect") abs(eff) else 0.8   # presence/direction fixed
  data.table(gene = g, modality = mod, direction = dir, disp_size = size)
}

rows <- list()
for (i in seq_len(nrow(core))) {
  r <- core[i]
  rows[[length(rows) + 1]] <- mk(r$human_symbol, "Human bulk",    r$lit_bulk,
                                 r$bulk_shrunk_logFC, "effect")
  rows[[length(rows) + 1]] <- mk(r$human_symbol, "Mouse bulk",    r$lit_mouse,
                                 r$mouse_meta_logFC, "effect")
  rows[[length(rows) + 1]] <- mk(r$human_symbol, "sc-Hepatocyte", r$lit_schep,
                                 r$sc_hepatocyte_logFC, "effect")
  rows[[length(rows) + 1]] <- mk(r$human_symbol, "Hotspot",       r$lit_hotspot,
                                 NA, "presence")
  rows[[length(rows) + 1]] <- mk(r$human_symbol, "Coloc PP.H4",   r$lit_coloc,
                                 NA, "presence")
}
long <- rbindlist(rows)
long[, modality := factor(modality, levels = mod_levels)]
long[, gene := factor(gene, levels = rev(core$human_symbol))]  # HKDC1 top

# ---------------------------------------------------------------------------
# Hero numbers
# ---------------------------------------------------------------------------
hk <- core[human_symbol == "HKDC1"]
hk_mods <- if (nrow(hk)) {
  mod_levels[c(hk$lit_bulk, hk$lit_mouse, hk$lit_schep, hk$lit_hotspot,
               hk$lit_coloc)]
} else character(0)
cat(sprintf("[hero] convergent core = %d genes (robust bulk DEG AND coloc PP.H4>0.5)\n",
            nrow(core)))
cat("[hero] genes:", paste(core$human_symbol, collapse = ", "), "\n")
cat(sprintf("[hero] HKDC1 lights %d modalities: %s (coloc PP.H4=%.3f)\n",
            ifelse(nrow(hk), hk$n_lit, 0L), paste(hk_mods, collapse = ", "),
            ifelse(nrow(hk), hk$coloc_best_pp4, NA_real_)))

fwrite(core, file.path(DATA_DIR, "crossmodal_convergence_matrix.csv"))

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
dir_colors <- c(Up = "#C9265E", Down = "#1565C0", Present = "#616161")

p <- ggplot(long, aes(x = modality, y = gene)) +
  geom_point(aes(fill = direction, size = disp_size),
             shape = 21, color = "white", stroke = 0.2) +
  scale_fill_manual(values = dir_colors, name = NULL,
                    breaks = c("Up", "Down", "Present")) +
  scale_size_continuous(range = c(1.4, 4.2), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Cross-modal convergent core") +
  theme_masld(base_size = 7) +
  theme(
    plot.title       = element_text(size = 7.3, face = "bold", margin = margin(b = 5)),
    axis.text.x.top  = element_text(size = 6, angle = 40, hjust = 0, face = "bold"),
    axis.text.y      = element_text(size = 6, face = "italic"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "bottom",
    legend.key.size  = unit(0.18, "cm")
  )

ggsave(OUT_PDF, p,
       width  = 78 / 25.4,
       height = max(70, 6 * nrow(core) + 28) / 25.4,
       units  = "in",
       device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
