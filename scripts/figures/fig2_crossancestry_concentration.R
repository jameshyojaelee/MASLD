#!/usr/bin/env Rscript
# fig2_crossancestry_concentration.R  (2026-07-05)
# Fig 2 (main) panel — "cross-ancestry probability-mass concentration".
#
# Illustrates how joint cross-ancestry fine-mapping redistributes posterior
# probability relative to an EUR-specific model. It is a model-resolution
# comparison, not identification of the causal variant.
#
# SOURCE (2026-07-05): rebuilt on the within-MVP N-way cross-ancestry run
#   (SuSiEx + meSuSiE; EUR/AFR/AMR/EAS; 7 MVP traits; 411 shared loci) that
#   REPLACES the retired enzyme-only EUR/EAS run. All PIPs come from the
#   authoritative _mvp gene summaries (susiex_max_pip / mesusie per-ancestry-combo
#   max_pip_*), NOT from combined_finemapping.csv, whose cross-ancestry
#   (susiex_pip / mesusie_pip) columns are stale/partial w.r.t. the MVP run.
#
# WHAT THIS PANEL SHOWS (and why it is the HONEST metric):
#   One horizontal dumbbell per established MASLD / liver locus (canonical set;
#   the full per-gene record is written to the source CSV).
#     left point  = EUR-specific lead PIP  = best PIP in an EUR-restricted
#                   meSuSiE credible set (max_pip_EUR) — the within-MVP EUR-only
#                   signal for that locus.
#     right point = joint lead PIP         = best PIP in a cross-ancestry SHARED
#                   meSuSiE credible set (max over the >=2-ancestry combo columns,
#                   e.g. max_pip_EUR_AFR_AMR).
#   These are DIFFERENT credible-set classes (neither a subset of the other), so
#   the EUR-specific -> shared delta can go UP or DOWN per locus — an honest test,
#   unlike a naive max(susiex, mesusie) which is >= the EUR-specific value by
#   construction. The SuSiEx overall joint PIP (susiex_max_pip) is retained as a
#   reference column in the source CSV. Arrow colored by direction (gain = teal,
#   loss/flat = grey); guides at PIP 0.5 and 0.9; rows ordered by joint PIP.
#
# DATA: GWAS/finemapping/results/susiex_mvp/susiex_gene_summary_mvp.csv
#       GWAS/finemapping/results/mesusie_mvp/mesusie_gene_summary_mvp.csv
#
# Output: figures/main/fig2_genetics/panels/FigS2J_crossancestry_pip_concentration.pdf
#         figures/main/fig2_genetics/panels/FigS2J_crossancestry_pip_concentration_source.csv
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- Sys.getenv("FIG2_SUPP_OUT_DIR", unset = file.path(FIG3_DIR, "panels"))
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Load the within-MVP N-way cross-ancestry gene summaries (authoritative).
#    combined_finemapping.csv is NOT used: its cross-ancestry columns are stale
#    w.r.t. the MVP run.
# ---------------------------------------------------------------------------
gs <- fread(file.path(BASE, "GWAS/finemapping/results/susiex_mvp/susiex_gene_summary_mvp.csv"),
            select = c("GeneSymbol", "susiex_max_pip", "locus_ids"))
mg <- fread(file.path(BASE, "GWAS/finemapping/results/mesusie_mvp/mesusie_gene_summary_mvp.csv"))
gs <- gs[GeneSymbol != "" & !is.na(GeneSymbol)]
mg <- mg[GeneSymbol != "" & !is.na(GeneSymbol)]

# meSuSiE per-ancestry-combo PIP columns. Split into single-ancestry-restricted
# (1 token, e.g. max_pip_EUR) vs cross-ancestry SHARED (>=2 tokens, e.g.
# max_pip_EUR_AFR_AMR). The honest concentration test compares the EUR-SPECIFIC
# credible-set PIP against the best SHARED credible-set PIP: these are DIFFERENT
# credible-set classes (neither a subset of the other), so the delta can go up
# OR down per locus -- unlike a naive max(susiex, mesusie) which is >= the
# EUR-specific value by construction (a tautological "gain" at every locus).
mp_cols  <- grep("^max_pip_", names(mg), value = TRUE)
mp_tok   <- vapply(sub("^max_pip_", "", mp_cols),
                   function(x) length(strsplit(x, "_")[[1]]), integer(1))
mp_multi <- mp_cols[mp_tok >= 2]                 # shared (>=2 ancestry) credible sets

# ---------------------------------------------------------------------------
# 2. One row per gene: EUR-specific baseline, best SHARED cross-ancestry PIP,
#    its ancestry composition, and the SuSiEx joint PIP (reference column).
# ---------------------------------------------------------------------------
res <- merge(gs, mg[, c("GeneSymbol", "mesusie_max_pip", mp_cols), with = FALSE],
             by = "GeneSymbol", all = TRUE)
res[, eur_pip      := suppressWarnings(as.numeric(max_pip_EUR))]       # EUR-specific CS PIP
res[, susiex_joint := suppressWarnings(as.numeric(susiex_max_pip))]    # SuSiEx overall joint (reference)

# best SHARED (cross-ancestry) credible-set PIP + which ancestry combo achieves it
shmat  <- as.matrix(res[, ..mp_multi]); shmat[is.na(shmat)] <- -Inf
sh_idx <- max.col(shmat, ties.method = "first")
res[, joint_pip := apply(shmat, 1, max)]
res[!is.finite(joint_pip), joint_pip := NA_real_]
res[, ancestries := gsub("_", "+", sub("^max_pip_", "", mp_multi[sh_idx]))]
res[is.na(joint_pip), ancestries := NA_character_]
res[, joint_method := "meSuSiE shared CS"]

# representative physical position (first locus_id tag; label/reference only)
first_lid <- sub(";.*", "", res$locus_ids)
mm <- regmatches(first_lid, regexec("chr([0-9XY]+)_([0-9]+)$", first_lid))
res[, chr := suppressWarnings(as.integer(sapply(mm, function(x) if (length(x) == 3) x[2] else NA)))]
res[, pos := suppressWarnings(as.integer(sapply(mm, function(x) if (length(x) == 3) x[3] else NA)))]
res[, gene := GeneSymbol]

# ---------------------------------------------------------------------------
# 3. Source CSV — full per-gene record with an EUR-specific baseline + joint PIP
# ---------------------------------------------------------------------------
res[, delta := joint_pip - eur_pip]
full <- res[!is.na(joint_pip) & !is.na(eur_pip)]
src <- full[, .(gene, chr, pos, ancestries, eur_pip = round(eur_pip, 4),
                joint_pip = round(joint_pip, 4), joint_method,
                susiex_joint_pip = round(susiex_joint, 4),
                delta = round(delta, 4))][order(-joint_pip)]
fwrite(src, file.path(PANEL_DIR, "FigS2J_crossancestry_pip_concentration_source.csv"))

# ---------------------------------------------------------------------------
# 4. Dumbbell plot — established MASLD / liver loci (the full set is in the CSV;
#    377 SuSiEx / 404 meSuSiE loci are too many to draw, so the panel highlights
#    the canonical loci, matching the panel's original density).
# ---------------------------------------------------------------------------
CANON <- c("GCKR", "PNPLA3", "TM6SF2", "MBOAT7", "HSD17B13", "TRIB1", "GPAM",
           "HNF1A", "HNF1B", "GATAD2A", "MLXIPL", "ABO", "HFE", "RECQL4",
           "GGT1", "SERPINA1", "MARC1", "APOE", "PPP1R3B", "TOR1B")
pl <- full[gene %in% CANON]
if (anyDuplicated(pl$gene)) pl <- pl[, .SD[which.max(joint_pip)], by = gene]
pl[, direction := fifelse(delta > 0.005, "gain", "loss/flat")]
setorder(pl, joint_pip, eur_pip)
# Use a NUMERIC row index (yi) on a continuous y-axis with manual gene labels.
# This lets the 0.5/0.9 header tics and the GCKR/GPAM callouts share the same
# coordinate system as the dumbbells (annotate() needs a numeric y).
pl[, yi := .I]

n_eur     <- nrow(pl)
n_gain    <- sum(pl$direction == "gain")
n_full    <- nrow(full)
gckr      <- pl[gene == "GCKR"]
gpam      <- pl[grepl("GPAM", gene)]

dir_cols <- c("gain" = "#00695C", "loss/flat" = "#9E9E9E")

p <- ggplot(pl, aes(y = yi)) +
  # 0.5 / 0.9 PIP guides
  geom_vline(xintercept = c(0.5, 0.9), linetype = "dashed",
             linewidth = 0.25, color = "grey75") +
  annotate("text", x = 0.5, y = n_eur + 0.85, label = "0.5", size = PUB_GEOM_TEXT,
           color = "black", vjust = 0) +
  annotate("text", x = 0.9, y = n_eur + 0.85, label = "0.9", size = PUB_GEOM_TEXT,
           color = "black", vjust = 0) +
  # connecting arrow EUR-only -> joint
  geom_segment(aes(x = eur_pip, xend = joint_pip, yend = yi, color = direction),
               linewidth = 0.55,
               arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  # EUR-only baseline point (open) and joint point (filled)
  geom_point(aes(x = eur_pip), shape = 21, fill = "white", color = "grey45",
             size = 1.6, stroke = 0.4) +
  geom_point(aes(x = joint_pip, color = direction), size = 1.7) +
  scale_color_manual(values = dir_cols, name = NULL,
                     labels = c(gain = "joint raised PIP",
                                "loss/flat" = "joint lowered / flat")) +
  scale_x_continuous(limits = c(0, 1.03), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(breaks = pl$yi, labels = pl$gene,
                     limits = c(0.4, n_eur + 1.5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Lead-variant PIP (EUR-specific ○ → joint cross-ancestry ●)",
    y = NULL
  ) +
  theme_masld() + theme_pub() +
  theme(
    # legend in the empty upper-left interior (top rows all sit at high PIP, so
    # the low-PIP side of those rows is blank)
    legend.position = c(0.02, 0.80),
    legend.justification = c(0, 0.5),
    legend.background = element_rect(fill = scales::alpha("white", 0.65), color = NA),
    legend.key.size = unit(0.18, "cm"),
    axis.text.y = element_text(face = "italic"),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(4, 8, 3, 3)
  )

# In-plot callouts for the two headline single-ancestry-weak gains (the loci
# that move farthest, to make the "concentration" message concrete). Placed
# above the joint point, away from the legend.
if (nrow(gpam) == 1)
  p <- p + annotate("text", x = gpam$joint_pip, y = gpam$yi,
                    label = sprintf("%.2f→%.2f", gpam$eur_pip, gpam$joint_pip),
                    size = PUB_GEOM_TEXT, color = "#00695C", vjust = -1.3, hjust = 0.6)
if (nrow(gckr) == 1)
  p <- p + annotate("text", x = gckr$eur_pip, y = gckr$yi,
                    label = sprintf("%.2f→%.2f", gckr$eur_pip, gckr$joint_pip),
                    size = PUB_GEOM_TEXT, color = "#00695C", vjust = -1.3, hjust = 1.08)

message(sprintf(
  "[caption] Joint fine-mapping redistributes probability mass. MVP cross-ancestry (EUR/AFR/AMR/EAS): joint inference raised the lead-variant PIP above the EUR-specific model at %d of %d displayed MASLD/liver loci; %d lowered, %d unchanged. This is method-resolution evidence, not proof of the causal variant or ancestry-specific biology.",
  n_gain, n_eur, sum(pl$delta < -0.005), sum(abs(pl$delta) <= 0.005)))

save_fig(p, file.path(PANEL_DIR, "FigS2J_crossancestry_pip_concentration.pdf"),
         width = fig_half_width, height = 2.9)

# ---------------------------------------------------------------------------
# 7. Console summary
# ---------------------------------------------------------------------------
cat("\n=== cross-ancestry PIP concentration (MVP N-way) ===\n")
cat(sprintf("Genes with a joint PIP + EUR-specific baseline (full source CSV): %d\n", n_full))
cat(sprintf("Plotted canonical MASLD/liver loci: %d; joint RAISED lead PIP above the EUR-specific signal in %d.\n",
            n_eur, n_gain))
cat("\nPer-locus (plotted canonical loci, ordered by joint PIP):\n")
print(pl[order(-joint_pip), .(gene, ancestries, eur_pip = round(eur_pip, 3),
                              joint_pip = round(joint_pip, 3), joint_method,
                              delta = round(delta, 3))])
cat("\nWrote:\n  ", file.path(PANEL_DIR, "FigS2J_crossancestry_pip_concentration.pdf"),
    "\n  ", file.path(PANEL_DIR, "FigS2J_crossancestry_pip_concentration_source.csv"), "\n")
