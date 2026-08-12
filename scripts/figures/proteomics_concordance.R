#!/usr/bin/env Rscript
# KEY MESSAGE: Atlas DEGs replicate at the protein level in TWO independent DIA-MS
# proteomes that sample different compartments — liver TISSUE (PXD051911) and
# PLASMA (PXD052937) — and in BOTH, mRNA-protein log2FC concordance RISES as the
# gene set is restricted from all measured genes -> primary DEGs -> the
# cross-species-conserved core. Conservation acts as a translatability filter:
#     liver   rho 0.34 -> 0.53 -> 0.61
#     plasma  rho 0.31 -> 0.50 -> 0.52
# (Koussounadis 2015 DEG-enrichment principle: mRNA-protein correlation is
# strongest among the genes that actually change.)
#
# DESIGN: paired per-gene scatters (Liver | Plasma), x = mRNA log2FC, y = protein
# log2FC; faded "all genes" cloud + DEGs (magenta) overlaid; Spearman rho + n
# annotated for each of those two gene-sets per compartment. NO lollipop
# (see memory/feedback-no-lollipop). rho/n recomputed from the per-gene table so
# the points and the annotated statistics are self-consistent; cross-checked
# against mrna_protein_concordance_stratified.csv.
#
# SIMPLIFIED 2026-07-02 (→ callout 4b): dropped the y=x dashed diagonal and the
# lm fit line (declutter), and removed the cross-species Conserved set from both
# the points and the annotation — conserved core has its own dedicated panel (4g),
# so here it only added a third indistinguishable teal cloud. Now two sets only:
# All genes vs Primary DEG.
#
# CAVEAT (legend): the two DIA-MS sources are NOT like-for-like replicates —
# liver tissue measures the biology directly; plasma is the non-invasive
# biomarker complement.
#
# Output: figures/main/fig4_validation/panels/figS4a.pdf (renamed 2026-07-07;
#         was fig4a_proteomics_concordance.pdf; also serves as Supp Fig 4a —
#         see docs/PAPER.md and archived Fig S4 rationale)
# Env:    rnaseq

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Per-gene mRNA(dream/C2)-protein concordance, both compartments ────────────
v3 <- fread(file.path(BASE,
  "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv"),
  select = c("gene", "protein_logFC", "bulk_logFC", "dataset", "dream_comparator"))
v3 <- v3[dream_comparator == "disease_vs_control" &
         dataset %in% c("PXD051911", "PXD052937") &
         is.finite(protein_logFC) & is.finite(bulk_logFC)]
v3[, compartment := fifelse(dataset == "PXD051911", "Liver tissue", "Plasma")]
v3[, compartment := factor(compartment, levels = c("Liver tissue", "Plasma"))]

# ── Gene-set membership: primary DEG (atlas C2) ───────────────────────────────
atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("human_symbol", "bulk_padj", "bulk_logFC"))
deg_genes <- atlas[!is.na(bulk_padj) & bulk_padj < 0.05 & abs(bulk_logFC) > 0.3, human_symbol]

# Category (DEG > other) drives colour + layering. (Conserved set removed 2026-07-02;
# conservation as a translatability filter has its own dedicated panel, 4g.)
v3[, category := fifelse(gene %in% deg_genes, "Primary DEG", "All genes")]
v3[, category := factor(category, levels = c("All genes", "Primary DEG"))]

# ── Spearman rho + n per compartment x gene-set (self-consistent w/ the points) ─
# "All genes" = the full overlap; DEG = the restricted subset.
rho_tab <- rbindlist(lapply(levels(v3$compartment), function(cp) {
  sub <- v3[compartment == cp]
  sets <- list("All genes"   = sub,
               "Primary DEG" = sub[gene %in% deg_genes])
  rbindlist(lapply(names(sets), function(s) {
    d <- sets[[s]]
    data.table(compartment = cp, set = s, n = nrow(d),
               rho = if (nrow(d) > 3) cor(d$bulk_logFC, d$protein_logFC, method = "spearman") else NA_real_)
  }))
}))
rho_tab[, compartment := factor(compartment, levels = c("Liver tissue", "Plasma"))]
rho_tab[, set := factor(set, levels = c("All genes", "Primary DEG"))]
# label text block per facet
lab_dt <- rho_tab[, .(label = paste(sprintf("%-9s ρ=%.2f (n=%s)",
                       set, rho, formatC(n, big.mark = ",", format = "d")),
                       collapse = "\n")), by = compartment]

# ── Colours + display clipping ───────────────────────────────────────────────
cat_cols <- c("All genes" = masld_colors$ns, "Primary DEG" = masld_colors$up)
LIM <- 4
v3[, `:=`(x = pmax(pmin(bulk_logFC, LIM), -LIM),
          y = pmax(pmin(protein_logFC, LIM), -LIM))]

# ── Plot: 2-facet scatter ────────────────────────────────────────────────────
p <- ggplot(v3, aes(x, y)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.25) +
  rasterize_layer(geom_point(data = v3[category == "All genes"],
                             color = cat_cols["All genes"], size = 0.18, alpha = 0.16, shape = 16), dpi = 600) +
  rasterize_layer(geom_point(data = v3[category == "Primary DEG"],
                             color = cat_cols["Primary DEG"], size = 0.3, alpha = 0.5, shape = 16), dpi = 600) +
  geom_text(data = lab_dt, aes(x = -LIM, y = LIM, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = GEOM_TEXT_6PT, family = "Helvetica",
            lineheight = 0.95, color = "black") +
  facet_wrap(~ compartment, nrow = 1) +
  coord_cartesian(xlim = c(-LIM, LIM), ylim = c(-LIM, LIM), clip = "off") +
  labs(x = expression("mRNA "*log[2]*"FC (disease vs control)"),
       y = expression("protein "*log[2]*"FC")) +
  # legend for the 3 gene-sets (manual, since colours are hard-set per layer)
  guides(color = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  theme_masld_compact() +
  theme(axis.text = element_text(color = "black"),
        panel.spacing.x = unit(0.3, "cm"))

# manual colour legend via dummy layer
p <- p + geom_point(data = data.frame(x = NA, y = NA, category = names(cat_cols)),
                    aes(color = category), na.rm = TRUE) +
  scale_color_manual(values = cat_cols, name = NULL,
                     breaks = names(cat_cols)) +
  theme(legend.position = "bottom", legend.key.size = unit(0.18, "cm"),
        legend.margin = margin(t = -5), plot.margin = margin(3, 4, 1, 3))

out <- file.path(FIG4_DIR, "panels", "figS4a.pdf")
save_fig(p, out, width = fig_full_width * 0.66, height = 1.9)
message("Saved: ", out)

# sidecar + cross-check vs canonical stratified table
fwrite(rho_tab, file.path(FIG4_DIR, "data", "proteomics_concordance.csv"))
message("[proteomics] rho/n per compartment x gene-set:")
print(rho_tab)
strat <- tryCatch(fread(file.path(BASE,
  "Analysis/Proteomics/results/mrna_protein_concordance_stratified.csv")), error = function(e) NULL)
if (!is.null(strat)) { message("[proteomics] canonical stratified table (cross-check):"); print(strat) }
