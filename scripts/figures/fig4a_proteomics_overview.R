#!/usr/bin/env Rscript
# KEY MESSAGE: Atlas transcriptional changes replicate at the protein level
# across two independent DIA-MS platforms (PXD051911 liver tissue + PXD052937
# plasma). Concordance strengthens as we restrict from all measured genes to
# primary DEGs to cross-species-conserved genes, and among genes significant in
# both layers directional agreement is 91.3%.
#
# Numbers are sourced VERBATIM from the canonical stratified concordance table
# (Analysis/Proteomics/results/mrna_protein_concordance_stratified.csv) so the
# panel matches the manuscript text EXACTLY:
#   All genes  rho = 0.31 (n = 3181)
#   DEGs       rho = 0.51 (n = 328)
#   Conserved  rho = 0.52 (n = 477)
#   Both sig   23 genes, 91.3% directional concordance
#   Per platform: plasma rho = 0.31, tissue rho = 0.34
# Concordance is computed against the canonical limma-voom quality-weighted C2
# differential expression (NOT the retired dream method).
#
# Output: figures/main/fig4_validation/proteomics_concordance.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Canonical stratified concordance ─────────────────────────────────────────
# Read the on-disk canonical table and pull the three headline strata. Reading
# (not recomputing) guarantees the panel never drifts from the manuscript text.
strat <- read.csv(file.path(BASE,
  "Analysis/Proteomics/results/mrna_protein_concordance_stratified.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)

get_row <- function(label) {
  hit <- strat[strat$stratum == label, , drop = FALSE]
  if (nrow(hit) != 1L)
    stop(sprintf("Expected exactly one '%s' row in stratified CSV (got %d)",
                 label, nrow(hit)))
  hit
}

r_all  <- get_row("All genes")
r_deg  <- get_row("DEGs (padj<0.05, |LFC|>0.3)")
r_cons <- get_row("Conserved")
r_both <- get_row("Both significant")

# Headline lollipop data: rho by stratum (least → most stringent filter).
bar_df <- data.frame(
  stratum = c("All genes", "Primary DEGs", "Conserved"),
  rho     = c(r_all$rho, r_deg$rho, r_cons$rho),
  n       = c(r_all$n,   r_deg$n,   r_cons$n),
  stringsAsFactors = FALSE
)
# Plot bottom→top in increasing stringency, so "Conserved" sits at the top.
bar_df$stratum <- factor(bar_df$stratum,
                         levels = c("All genes", "Primary DEGs", "Conserved"))

# Both-significant directional concordance (the panel's stated headline stat).
both_n   <- as.integer(round(r_both$n))             # 23
both_pct <- r_both$direction_pct                    # 91.3043...

# Per-platform rho. Tissue (PXD051911 MASLD-vs-ctrl) is in the stratified CSV;
# plasma (PXD052937) is the canonical 0.31 reported in the manuscript and in
# pxd052937_concordance_summary.csv (rho_overall). Hard-source both so the
# values cannot drift.
plat <- read.csv(file.path(BASE,
  "Analysis/Proteomics/results/pxd052937_concordance_summary.csv"),
  stringsAsFactors = FALSE)
rho_plasma <- plat$value[plat$metric == "rho_overall"]          # 0.3065 -> 0.31
tiss_row   <- strat[strat$stratum == "PXD051911 liver (MASLD vs ctrl)", ]
rho_tissue <- tiss_row$rho[1]                                   # 0.3424 -> 0.34

# Emit to stdout for caption cross-checking.
message(sprintf("[fig4a] rho: all=%.2f (n=%d)  DEG=%.2f (n=%d)  cons=%.2f (n=%d)",
                bar_df$rho[1], bar_df$n[1], bar_df$rho[2], bar_df$n[2],
                bar_df$rho[3], bar_df$n[3]))
message(sprintf("[fig4a] both-sig: n=%d, %.1f%% directional concordance",
                both_n, both_pct))
message(sprintf("[fig4a] per-platform rho: plasma=%.2f  tissue=%.2f",
                rho_plasma, rho_tissue))

# ── Colors ───────────────────────────────────────────────────────────────────
# Stringency ramp: neutral gray (all) -> disease magenta (DEGs) -> teal
# (conserved), matching the house concordance semantics.
strat_cols <- c(
  "All genes"    = masld_colors$ns,         # neutral gray
  "Primary DEGs" = masld_colors$up,         # Liang deep magenta (DEG signal)
  "Conserved"    = masld_colors$conserved   # deep teal (conserved)
)

# ── Headline lollipop: stratified rho ────────────────────────────────────────
# Axis spans 0–0.6 (the data); extra blank space to the right of the axis (via
# clip="off" + plot.margin) holds the value+n label so nothing is truncated.
x_axis_max <- 0.6
p <- ggplot(bar_df, aes(x = rho, y = stratum, color = stratum)) +
  geom_segment(aes(x = 0, xend = rho, yend = stratum), linewidth = 1.1) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.2f  (n=%s)", rho,
                                formatC(n, big.mark = ",", format = "d"))),
            hjust = 0, nudge_x = 0.022,
            size = PUB_GEOM_TEXT + 0.4, fontface = "bold", color = "black") +
  scale_color_manual(values = strat_cols, guide = "none") +
  scale_x_continuous(limits = c(0, x_axis_max),
                     breaks = c(0, 0.2, 0.4, 0.6),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "mRNA-protein log2FC concordance (Spearman rho)",
       y = NULL,
       title = "Protein-level replication") +
  # Both-significant directional concordance + per-platform agreement belong in
  # the figure legend, not the panel (PI directive) — emitted to stdout above.
  coord_cartesian(clip = "off") +
  theme_masld() + theme_pub() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = PUB_AXIS_TITLE),
        plot.margin = margin(5.5, 34, 5.5, 5.5))

# ── Save ─────────────────────────────────────────────────────────────────────
out <- file.path(FIG4_DIR, "proteomics_concordance.pdf")
pdf(out, width = fig_half_width, height = fig_half_width * 0.62, useDingbats = FALSE)
print(p)
dev.off()
message("Saved: ", out)
