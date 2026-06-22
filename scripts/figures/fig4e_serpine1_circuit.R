#!/usr/bin/env Rscript
# KEY MESSAGE: Among fibroblast-derived secreted ligands, SERPINE1 stands out by
# pairing a high LIANA ligand score with strong bulk-disease upregulation — and
# its spatial structure strengthens sharply in steatotic liver. A blood-accessible
# fibroblast signal.
#
# NOTE (provenance): SERPINE1's LIANA fibroblast-ligand score is 0.916. It is NOT
# the single highest fibroblast-ligand score (LAMA4 0.972, NUCB2 0.968 rank above
# it); what distinguishes SERPINE1 is the COMBINATION of high ligand score AND
# high bulk LFC (+0.54) — LAMA4/NUCB2 are near-flat in bulk. The panel shows that
# combination honestly rather than asserting rank #1.
#
# Panels:
#   i.  Fibroblast-derived ligands: bulk LFC (x) x LIANA score (y), SERPINE1 hero
#   ii. SERPINE1 spatial autocorrelation, healthy -> steatotic (GSE192741)
#
# Bulk numbers are the canonical limma-voom quality-weighted C2 (bulk_*), NOT the
# retired dream method. Stats live in the figure legend (emitted to stdout).
#
# Output: figures/main/fig4_validation/serpine1_ligand.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data ──────────────────────────────────────────────────────────────────────
secretome <- read.csv(file.path(BASE,
  "RNA-seq/results/secretome_chain/secretome_chain_triple_concordant_MASLD_up.csv"),
  stringsAsFactors = FALSE) %>%
  mutate(gene = ligand_complex)

# Per-condition spatial autocorrelation (single source of truth; GSE192741).
autocorr_h <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Healthy.csv"),
  stringsAsFactors = FALSE)
autocorr_s <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Steatotic.csv"),
  stringsAsFactors = FALSE)
serp_moran_h <- autocorr_h$C[autocorr_h$Gene == "SERPINE1"][1]   # 0.059 healthy
serp_moran_s <- autocorr_s$C[autocorr_s$Gene == "SERPINE1"][1]   # 0.541 steatotic

# ── Panel i: fibroblast-derived ligand landscape ─────────────────────────────
# Restrict to fibroblast-sourced ligands (the manuscript's scope). One row per
# gene (best LIANA score). SERPINE1 is the only highlighted point.
fib <- secretome %>%
  mutate(src = ifelse(!is.na(primary_source_ct) & primary_source_ct != "",
                      primary_source_ct, a1_primary)) %>%
  filter(grepl("fibro", src, ignore.case = TRUE),
         !is.na(max_score_diff), !is.na(bulk_lfc)) %>%
  group_by(gene) %>% slice_max(max_score_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(is_focal = gene == "SERPINE1")

p_land <- ggplot(fib, aes(x = bulk_lfc, y = max_score_diff)) +
  geom_point(aes(color = is_focal, size = is_focal)) +
  geom_text_repel(aes(label = gene, fontface = ifelse(is_focal, "bold.italic", "italic"),
                      color = is_focal),
                  size = PUB_GEOM_TEXT, box.padding = 0.3, segment.size = 0.2,
                  max.overlaps = Inf, seed = 1, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = masld_colors[["up"]], `FALSE` = "#C9A0B4"),
                     guide = "none") +
  scale_size_manual(values = c(`TRUE` = 3.2, `FALSE` = 1.4), guide = "none") +
  scale_y_continuous(limits = c(0.65, 1.0)) +
  labs(x = "Bulk mRNA log2FC", y = "LIANA ligand score",
       title = "Fibroblast-derived ligands") +
  theme_masld() + theme_pub()

# ── Panel ii: SERPINE1 spatial autocorrelation, healthy -> steatotic ─────────
spat_df <- data.frame(
  cond  = factor(c("Healthy", "Steatotic"), levels = c("Healthy", "Steatotic")),
  moran = c(serp_moran_h, serp_moran_s),
  fill  = c(masld_colors[["control"]], masld_colors[["up"]])
)
p_spat <- ggplot(spat_df, aes(x = cond, y = moran, fill = fill)) +
  geom_col(width = 0.6) +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, max(spat_df$moran) * 1.08),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Spatial autocorrelation\n(Moran's I)",
       title = "SERPINE1 spatial") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(size = PUB_AXIS_TEXT + 0.5))

# ── Assemble ──────────────────────────────────────────────────────────────────
serp <- fib %>% filter(is_focal)
message(sprintf("[serpine1 legend] SERPINE1 fibroblast LIANA score = %.3f; bulk logFC = %.3f; Moran's I %.3f (healthy) -> %.3f (steatotic)",
                serp$max_score_diff[1], serp$bulk_lfc[1], serp_moran_h, serp_moran_s))

p_out <- (p_land | p_spat) + plot_layout(widths = c(2.2, 1))

out <- file.path(FIG4_DIR, "serpine1_ligand.pdf")
pdf(out, width = fig_full_width * 0.92, height = 2.5, useDingbats = FALSE)
print(p_out)
dev.off()
message("Saved: ", out)
