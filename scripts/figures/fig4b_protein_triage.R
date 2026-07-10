#!/usr/bin/env Rscript
# ==============================================================================
# Fig 4b (HEADLINE) — protein-layer triage of the prioritized MASLD target set.
# Across the WHOLE prioritized universe measured in an independent liver proteomics
# cohort (PXD051911 liver DIA-MS, n=58), mRNA rank only weakly predicts protein rank.
# This panel shows the full mRNA-vs-protein log2FC distribution and partitions the
# measured prioritized set into 3 classes, handing forward the protein-confirmed core.
#   scatter: x = mRNA log2FC (pooled bulk RNA-seq), y = protein log2FC (PXD051911)
#            colour = 3-way class (prioritized) vs light-gray rest-of-proteome context
#   bar    : the prioritized measured set partitioned into the 3 classes (counts + %)
# NOT a pass/fail filter: the majority "not yet protein-corroborated" class is buffering
# AND/OR below DIA-MS detection depth, NOT evidence the target is wrong.
# Output: figures/main/fig4_validation/panels/fig4b_protein_triage.pdf
# Env: rnaseq
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
set.seed(42)
FAM  <- "Helvetica"
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── data: liver DIA-MS mRNA↔protein concordance, prioritized universe ─────────
con <- fread(file.path(BASE, "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv"))
con <- unique(con[dataset == "PXD051911" & is.finite(bulk_logFC) & is.finite(protein_logFC)], by = "gene")
uni <- trimws(readLines(file.path(BASE, "Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt")))
uni <- uni[uni != ""]
con[, in_uni := gene %in% uni]

# whole-proteome vs prioritized-universe correlation (Spearman) + direction concordance
sgn_conc <- function(a, b) mean(sign(a) == sign(b))
rho_full  <- cor(con$bulk_logFC, con$protein_logFC, method = "spearman")
conc_full <- sgn_conc(con$bulk_logFC, con$protein_logFC)
du <- con[in_uni == TRUE]
rho_uni   <- cor(du$bulk_logFC, du$protein_logFC, method = "spearman")
conc_uni  <- sgn_conc(du$bulk_logFC, du$protein_logFC)

# ── 3-way classification of the PRIORITIZED measured set ──────────────────────
con[, cls := fifelse(!in_uni, "rest",
              fifelse(sign(bulk_logFC) == sign(protein_logFC) & protein_padj < 0.05, "confirmed",
              fifelse(sign(bulk_logFC) != sign(protein_logFC) & protein_padj < 0.05, "discordant",
              "not_corr")))]
n_uni   <- nrow(du)
n_conf  <- sum(con$cls == "confirmed")
n_disc  <- sum(con$cls == "discordant")
n_notc  <- sum(con$cls == "not_corr")
n_rest  <- sum(con$cls == "rest")
cat(sprintf("[fig4b-triage] measured=%d (universe=%d, rest=%d) | full rho=%.3f conc=%.1f%% | uni rho=%.3f conc=%.1f%%\n",
            nrow(con), n_uni, n_rest, rho_full, 100*conc_full, rho_uni, 100*conc_uni))
cat(sprintf("[fig4b-triage] 3-way (universe): confirmed=%d (%.1f%%) not_corroborated=%d (%.1f%%) discordant=%d (%.1f%%)\n",
            n_conf, 100*n_conf/n_uni, n_notc, 100*n_notc/n_uni, n_disc, 100*n_disc/n_uni))

# colours / draw order / point styling per class
COL <- c(rest = "#E0E0E0", not_corr = "#9E9E9E", confirmed = "#00695C", discordant = "#E91E63")
SZ  <- c(rest = 0.45,      not_corr = 0.7,       confirmed = 1.15,      discordant = 1.15)
AL  <- c(rest = 0.35,      not_corr = 0.55,      confirmed = 0.95,      discordant = 0.95)
lab_map <- c(rest = "Background proteome", not_corr = "Not significant",
             confirmed = "Significant, concordant", discordant = "Discordant")
con[, cls := factor(cls, levels = c("rest", "not_corr", "confirmed", "discordant"))]
setorder(con, cls)   # draw rest -> not_corr -> confirmed -> discordant (highlights on top)

# ── scatter ───────────────────────────────────────────────────────────────────
# single essential data annotation: the headline correlation for the prioritized set
# (full-proteome rho and cross-cohort/buffering narrative live in the message() CAPTION).
ann <- sprintf("ρ = %.2f · %.0f%% direction-concordant", rho_uni, 100*conc_uni)
xr <- range(con$bulk_logFC); yr <- range(con$protein_logFC)
sc <- ggplot(con, aes(bulk_logFC, protein_logFC, colour = cls, size = cls, alpha = cls)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.22, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.22, colour = "grey70") +
  rasterize_layer(geom_point(stroke = 0), dpi = 400) +
  # OLS best-fit over the prioritized universe (matches the reported rho_uni annotation)
  geom_smooth(data = du, aes(bulk_logFC, protein_logFC), method = "lm", formula = y ~ x,
              se = FALSE, linetype = "dashed", linewidth = 0.4, colour = house_ink,
              inherit.aes = FALSE) +
  scale_colour_manual(values = COL, labels = lab_map, name = NULL,
                      guide = guide_legend(override.aes = list(size = 1.4, alpha = 1), ncol = 1)) +
  scale_size_manual(values = SZ, guide = "none") +
  scale_alpha_manual(values = AL, guide = "none") +
  annotate("text", x = xr[1], y = yr[2], hjust = 0, vjust = 1, size = 6/.pt, family = FAM,
           label = ann) +
  labs(x = expression("mRNA "*log[2]*"FC (bulk RNA-seq)"),
       y = expression("protein "*log[2]*"FC (PXD051911 liver DIA-MS)")) +
  theme_masld() +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.size = unit(0.28, "lines"), legend.text = element_text(size = 6),
        legend.margin = margin(1, 2, 1, 2), plot.margin = margin(2, 1, 1, 2))

# ── 3-way VERTICAL stacked bar (prioritized measured set) — NOT a pass/fail filter ─────
# bottom→top: confirmed, discordant, not_corr. Every segment labelled to its right,
# so counts+% sit on their own line (no cramming — the vertical axis gives room).
bar <- data.table(
  cls = factor(c("confirmed", "discordant", "not_corr"),
               levels = c("confirmed", "discordant", "not_corr")),
  n   = c(n_conf, n_disc, n_notc))
setorder(bar, cls)
bar[, ymax := cumsum(n)][, ymin := ymax - n][, mid := (ymin + ymax)/2]
bar[, pct := 100 * n / sum(n)]
# data labels only (class name + count/%); mechanism narrative lives in the CAPTION
bar_lab <- c(confirmed = sprintf("Significant,\nconcordant\n%d (%.1f%%)", n_conf, 100*n_conf/n_uni),
             discordant = sprintf("Discordant\n%d (%.1f%%)", n_disc, 100*n_disc/n_uni),
             not_corr  = sprintf("Not significant\n%s (%.0f%%)",
                                 format(n_notc, big.mark = ","), 100*n_notc/n_uni))
# the two tiny segments sit near y=0 → spread their labels out with leader lines;
# the big majority labels at its own mid.
lab_x  <- 1.22                    # label column (bar occupies x = 0..1)
lead_x <- 1.14                    # leader elbow
y_conf <- n_uni * 0.10            # spread-out target y for confirmed label
y_disc <- n_uni * 0.30            # spread-out target y for discordant label
bp <- ggplot(bar) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = cls),
            colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = COL[c("confirmed","discordant","not_corr")], guide = "none") +
  # big neutral majority: label to the right at its own mid
  annotate("text", x = lab_x, y = bar[cls=="not_corr", mid], label = bar_lab["not_corr"],
           size = 6/.pt, family = FAM, colour = house_ink, hjust = 0, vjust = 0.5, lineheight = 0.9) +
  # tiny confirmed: leader from segment mid out to spread-out label
  annotate("segment", x = 1, xend = lead_x, y = bar[cls=="confirmed", mid], yend = y_conf,
           linewidth = 0.25, colour = "#00695C") +
  annotate("text", x = lab_x, y = y_conf, label = bar_lab["confirmed"],
           size = 6/.pt, family = FAM, colour = house_ink, hjust = 0, vjust = 0.5, lineheight = 0.9) +
  # tiny discordant: leader from segment mid out to spread-out label
  annotate("segment", x = 1, xend = lead_x, y = bar[cls=="discordant", mid], yend = y_disc,
           linewidth = 0.25, colour = "#E91E63") +
  annotate("text", x = lab_x, y = y_disc, label = bar_lab["discordant"],
           size = 6/.pt, family = FAM, colour = house_ink, hjust = 0, vjust = 0.5, lineheight = 0.9) +
  scale_x_continuous(limits = c(0, 2.55), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, n_uni), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(y = sprintf("prioritized targets measured\nin liver proteome (n = %s)", format(n_uni, big.mark = ",")),
       x = NULL) +
  theme_void() +
  theme(text = element_text(family = FAM, face = "plain"),
        axis.title.y = element_text(size = 6, family = FAM, angle = 90, margin = margin(r = 2)),
        plot.margin = margin(t = 2, r = 36, b = 1, l = 1))

# ── compose (scatter left, vertical bar right) ─────────────────────────────────
p <- (sc | bp) + plot_layout(widths = c(3.5, 1.35))
out <- file.path(FIG4_DIR, "panels", "fig4b_protein_triage.pdf")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
ggsave(out, p, width = 3.56, height = 2.13, device = grDevices::cairo_pdf)
cat("[fig4b-triage] saved:", out, "\n")

message(sprintf(paste0(
  "CAPTION (Fig 4b, HEADLINE): Protein-layer triage of the prioritized MASLD target set in an INDEPENDENT liver ",
  "proteomics cohort (PXD051911 liver DIA-MS, n=58 patients — a DIFFERENT cohort from the bulk RNA-seq, so the ",
  "mRNA-vs-protein log2FC correlation is attenuated by cross-cohort measurement noise ON TOP OF any real ",
  "post-transcriptional regulation). Each point is a gene measured at both layers; x = mRNA log2FC (pooled bulk ",
  "RNA-seq), y = protein log2FC (PXD051911). Across the full measured proteome (n=%s) Spearman rho=%.2f and %.0f%% ",
  "of genes agree in direction; restricted to the prioritized universe (n=%s measured of the %s-gene shared Fig 4 ",
  "target set = transcriptomic 8,088 union genetic 3,038) rho=%.2f and %.0f%% agree ",
  "-- both WELL ABOVE the 50%% chance level, so real shared mRNA->protein signal exists even though it is weak, ",
  "but mRNA rank only partially predicts protein rank. The prioritized measured set partitions into: ",
  "protein-confirmed = direction-concordant AND protein FDR<0.05 (n=%d, %.1f%%, handed forward as the ",
  "protein-corroborated core); discordant = opposite direction AND protein-significant (n=%d, %.1f%%); and ",
  "'not yet protein-corroborated' = the majority (n=%s, %.0f%%) -- buffering AND/OR below DIA-MS detection depth, ",
  "NOT evidence the target is wrong (this is a triage, not a pass/fail filter). Fig 4c independently shows these ",
  "liver protein measurements track histological severity, confirming the protein layer is biologically real."),
  format(nrow(con), big.mark=","), rho_full, 100*conc_full,
  format(n_uni, big.mark=","), format(length(uni), big.mark=","), rho_uni, 100*conc_uni,
  n_conf, 100*n_conf/n_uni, n_disc, 100*n_disc/n_uni, format(n_notc, big.mark=","), 100*n_notc/n_uni))
