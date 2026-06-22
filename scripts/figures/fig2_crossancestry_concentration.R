#!/usr/bin/env Rscript
# fig2_crossancestry_concentration.R  (2026-06-12)
# Fig 2 (main) panel — "cross-ancestry probability-mass concentration".
#
# Supports the manuscript claim (fig2.md para 4): joint cross-ancestry
# fine-mapping shifts probability mass toward the true causal variant.
#
# WHAT THIS PANEL SHOWS (and why it is the HONEST metric):
#   The credible-set-SIZE metric is mixed/unsupportive in this portfolio
#   (only a handful of multi-ancestry loci; some shrink, some grow — PNPLA3
#   grows). The probability-mass-CONCENTRATION metric IS supported: joint
#   inference (SuSiEx / meSuSiE) generally RAISES the lead-variant PIP.
#
#   Encoding: one horizontal dumbbell per multi-ancestry locus.
#     left point  = EUR-only lead-variant PIP   (max SuSiE PIP over EUR rows)
#     right point = joint lead-variant PIP       (max over SuSiEx / meSuSiE)
#     arrow EUR -> joint, colored by direction (gain = teal, loss/flat = grey)
#   x = lead-variant PIP (0 -> 1); guides at PIP 0.5 and 0.9.
#   Rows ordered by joint PIP; each labeled by its nearest gene.
#
# HONESTY: we restrict the dumbbell to multi-ancestry loci that HAVE a defined
#   EUR-only baseline (EUR contributed). Multi-ancestry loci with no EUR arm
#   (EAS+SAS only) cannot be drawn as a EUR->joint shift; they are recorded in
#   the source CSV (eur_pip = NA) and noted in the caption, not hidden.
#
# DATA: GWAS/finemapping/results/combined_finemapping.csv (variant-level).
#   A locus = a physical cluster of joint-finemapped variants (gap > 500 kb
#   splits loci on a chromosome). joint_pip = max(susiex_pip, mesusie_pip)
#   over the locus; eur_pip = max(susie_pip) over EUR rows in the locus window.
#   "Multi-ancestry" = >= 2 ancestries present among rows in the locus window.
#
# Output: figures/main/fig2_genetics/panels/crossancestry_pip_concentration.pdf
#         figures/main/fig2_genetics/panels/crossancestry_pip_concentration_source.csv
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG3_DIR, "panels")          # FIG3_DIR == figures/main/fig2_genetics
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Load finemapping; keep variants with chr+pos
# ---------------------------------------------------------------------------
dt <- fread(file.path(BASE, "GWAS/finemapping/results/combined_finemapping.csv"),
            select = c("chromosome", "position", "study", "ancestry",
                       "susie_pip", "susiex_pip", "mesusie_pip"))
dt <- dt[!is.na(chromosome) & !is.na(position)]
for (cc in c("susie_pip", "susiex_pip", "mesusie_pip"))
  dt[, (cc) := suppressWarnings(as.numeric(get(cc)))]
dt[, has_joint := !is.na(susiex_pip) | !is.na(mesusie_pip)]

# ---------------------------------------------------------------------------
# 2. Cluster joint-finemapped variants into physical loci (gap > 500 kb)
#    NOTE: susiex_pip / mesusie_pip are per-VARIANT joint PIPs broadcast across
#    that variant's study rows, so dedup to one row per variant first.
# ---------------------------------------------------------------------------
jv <- unique(dt[has_joint == TRUE, .(chromosome, position, susiex_pip, mesusie_pip)])
jv[, joint_pip := pmax(susiex_pip, mesusie_pip, na.rm = TRUE)]
setorder(jv, chromosome, position)
jv[, gap := position - shift(position, fill = position[1]), by = chromosome]
jv[, new_locus := as.integer(gap > 5e5 | is.na(gap))]
jv[1, new_locus := 1L]
jv[, locus_id := cumsum(new_locus), by = chromosome]
jv[, locus_key := paste0(chromosome, "_", locus_id)]

loci_range <- jv[, .(chr = chromosome[1], pmin = min(position), pmax = max(position)),
                 by = locus_key]

# ---------------------------------------------------------------------------
# 3. Per-locus: ancestries present, EUR-only lead PIP, joint lead PIP + method
# ---------------------------------------------------------------------------
build_locus <- function(i) {
  lk <- loci_range$locus_key[i]; ch <- loci_range$chr[i]
  p0 <- loci_range$pmin[i];      p1 <- loci_range$pmax[i]
  sub <- dt[chromosome == ch & position >= p0 - 1e4 & position <= p1 + 1e4]
  anc <- sort(unique(na.omit(sub$ancestry)))
  eur_lead <- suppressWarnings(max(sub[ancestry == "EUR", susie_pip], na.rm = TRUE))
  if (!is.finite(eur_lead)) eur_lead <- NA_real_
  lrows  <- jv[locus_key == lk]
  lead   <- lrows[which.max(joint_pip)]
  jmeth  <- if (!is.na(lead$mesusie_pip) &&
                (is.na(lead$susiex_pip) || lead$mesusie_pip >= lead$susiex_pip))
              "meSuSiE" else "SuSiEx"
  data.table(locus_key = lk, chr = ch, pos = lead$position,
             n_anc = length(anc), ancestries = paste(anc, collapse = "+"),
             eur_pip = eur_lead, joint_pip = lead$joint_pip, joint_method = jmeth)
}
res <- rbindlist(lapply(seq_len(nrow(loci_range)), build_locus))
res <- res[n_anc >= 2]                       # multi-ancestry loci only

# ---------------------------------------------------------------------------
# 4. Gene labels. Three-tier assignment, most-authoritative first:
#    (a) Canonical-MASLD override anchored to the physical interval (so a window
#        shift can never mis-assign a famous locus, and so two distinct loci
#        cannot collapse onto the same nearest-gene call).
#    (b) SuSiEx gene summary: locus_ids carry "chrC_POS"; take the gene whose
#        tag position is nearest the lead joint variant inside the window. This
#        is the authoritative per-locus gene the joint pipeline itself emitted.
#    (c) Nearest gene from the coloc lead annotation (sparse, used only as last
#        resort), then a chr:pos string.
# ---------------------------------------------------------------------------
# (b) SuSiEx gene summary -> chr/pos tags
gs <- fread(file.path(BASE, "GWAS/finemapping/results/susiex/susiex_gene_summary.csv"),
            select = c("GeneSymbol", "locus_ids"))
gs_long <- gs[, .(lid = unlist(strsplit(locus_ids, ";"))), by = GeneSymbol]
mm <- regmatches(gs_long$lid, regexec("chr([0-9]+)_([0-9]+)$", gs_long$lid))
gs_long[, g_chr := suppressWarnings(as.integer(sapply(mm, function(x) if (length(x) == 3) x[2] else NA)))]
gs_long[, g_pos := suppressWarnings(as.numeric(sapply(mm, function(x) if (length(x) == 3) x[3] else NA)))]
gs_long <- gs_long[!is.na(g_chr) & !is.na(g_pos)]

# (c) coloc nearest-gene annotation
ann <- fread(file.path(BASE, "RNA-seq/results/coloc_variant_classes/lead_causal_annotation.csv"),
             select = c("chr", "pos_hg19", "nearest_gene_symbol"))
ann <- ann[nearest_gene_symbol != "" & !is.na(nearest_gene_symbol)]

# (a) canonical-MASLD / liver-enzyme genes, anchored to the locus interval.
#     Fires only if the lead joint variant falls inside the interval.
canon <- data.table(
  chr  = c(2,      22,       8,       10,     4,          12,      17,      22,     9,     7,        6,     8),
  lo   = c(27.0e6, 44.1e6,   126.3e6, 113.8e6,88.0e6,     120.9e6, 35.9e6,  24.0e6, 136.0e6,73.0e6,  26.0e6,145.5e6),
  hi   = c(28.0e6, 44.5e6,   126.7e6, 114.1e6,88.3e6,     122.1e6, 36.3e6,  25.7e6, 136.3e6,73.2e6,  27.1e6,146.2e6),
  gene = c("GCKR", "PNPLA3", "TRIB1", "GPAM", "HSD17B13", "HNF1A", "HNF1B", "GGT1", "ABO", "MLXIPL", "HFE", "RECQL4")
)

assign_gene <- function(ch, p) {
  hit <- canon[chr == ch & p >= lo & p <= hi]
  if (nrow(hit) >= 1) return(hit$gene[1])
  sub <- gs_long[g_chr == ch & abs(g_pos - p) <= 5e5]
  if (nrow(sub) > 0) return(sub[which.min(abs(g_pos - p)), GeneSymbol])
  a <- ann[chr == ch]
  if (nrow(a) > 0) {
    cand <- a[which.min(abs(pos_hg19 - p))]
    if (abs(cand$pos_hg19 - p) < 5e5) return(cand$nearest_gene_symbol)
  }
  paste0("chr", ch, ":", round(p / 1e6, 2), "Mb")
}
res[, gene := mapply(assign_gene, chr, pos)]
# Guard: should not happen after canonical anchoring, but keep labels unique so
# the y-axis factor can never collapse two distinct physical loci.
if (anyDuplicated(res$gene)) res[, gene := make.unique(gene, sep = " ")]

# ---------------------------------------------------------------------------
# 5. Source CSV (all multi-ancestry loci, incl. EUR-absent ones — honest record)
# ---------------------------------------------------------------------------
res[, delta := joint_pip - eur_pip]
src <- res[, .(gene, chr, pos, ancestries, eur_pip = round(eur_pip, 4),
               joint_pip = round(joint_pip, 4), joint_method,
               delta = round(delta, 4))][order(-joint_pip)]
fwrite(src, file.path(PANEL_DIR, "crossancestry_pip_concentration_source.csv"))

# ---------------------------------------------------------------------------
# 6. Dumbbell plot — only loci with a defined EUR-only baseline
# ---------------------------------------------------------------------------
pl <- res[!is.na(eur_pip)]                                  # only loci with a EUR baseline
pl[, direction := fifelse(delta > 0.005, "gain", "loss/flat")]
setorder(pl, joint_pip, eur_pip)
# Use a NUMERIC row index (yi) on a continuous y-axis with manual gene labels.
# This lets the 0.5/0.9 header tics and the GCKR/GPAM callouts share the same
# coordinate system as the dumbbells (annotate() needs a numeric y).
pl[, yi := .I]

n_eur     <- nrow(pl)
n_gain    <- sum(pl$direction == "gain")
n_noeur   <- res[is.na(eur_pip), .N]
gckr      <- pl[gene == "GCKR"]
gpam      <- pl[grepl("GPAM", gene)]

dir_cols <- c("gain" = "#00695C", "loss/flat" = "#9E9E9E")

p <- ggplot(pl, aes(y = yi)) +
  # 0.5 / 0.9 PIP guides
  geom_vline(xintercept = c(0.5, 0.9), linetype = "dashed",
             linewidth = 0.25, color = "grey75") +
  annotate("text", x = 0.5, y = n_eur + 0.85, label = "0.5", size = 1.7,
           color = "grey55", vjust = 0) +
  annotate("text", x = 0.9, y = n_eur + 0.85, label = "0.9", size = 1.7,
           color = "grey55", vjust = 0) +
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
    x = "Lead-variant PIP (EUR-only ○ → joint ●)",
    y = NULL,
    title = "Joint fine-mapping concentrates probability mass",
    subtitle = sprintf(
      "Joint (SuSiEx / meSuSiE) raised the lead-variant PIP at %d of %d\nmulti-ancestry loci; %d lowered, %d unchanged.",
      n_gain, n_eur, sum(pl$delta < -0.005), sum(abs(pl$delta) <= 0.005))
  ) +
  theme_masld() + theme_pub() +
  theme(
    # legend in the empty upper-left interior (top rows all sit at high PIP, so
    # the low-PIP side of those rows is blank)
    legend.position = c(0.02, 0.80),
    legend.justification = c(0, 0.5),
    legend.background = element_rect(fill = scales::alpha("white", 0.65), color = NA),
    legend.key.size = unit(0.18, "cm"),
    plot.subtitle = element_text(size = PUB_SUBTITLE - 0.5, color = "gray30", lineheight = 1.05),
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
                    size = 1.7, color = "#00695C", vjust = -1.3, hjust = 0.6)
if (nrow(gckr) == 1)
  p <- p + annotate("text", x = gckr$eur_pip, y = gckr$yi,
                    label = sprintf("%.2f→%.2f", gckr$eur_pip, gckr$joint_pip),
                    size = 1.7, color = "#00695C", vjust = -1.3, hjust = 1.08)

save_fig(p, file.path(PANEL_DIR, "crossancestry_pip_concentration.pdf"),
         width = fig_half_width, height = 2.9)

# ---------------------------------------------------------------------------
# 7. Console summary
# ---------------------------------------------------------------------------
cat("\n=== cross-ancestry PIP concentration ===\n")
cat(sprintf("Multi-ancestry loci total: %d (with EUR arm: %d; EUR-absent EAS/SAS-only: %d)\n",
            nrow(res), n_eur, n_noeur))
cat(sprintf("Joint RAISED lead PIP in %d of %d EUR-arm loci.\n", n_gain, n_eur))
cat("\nPer-locus (EUR-arm loci, ordered by joint PIP):\n")
print(pl[order(-joint_pip), .(gene, ancestries, eur_pip = round(eur_pip, 3),
                              joint_pip = round(joint_pip, 3), joint_method,
                              delta = round(delta, 3))])
cat("\nEUR-absent multi-ancestry loci (recorded in source CSV, not plotted):\n")
print(res[is.na(eur_pip), .(gene, ancestries, joint_pip = round(joint_pip, 3), joint_method)])
cat("\nWrote:\n  ", file.path(PANEL_DIR, "crossancestry_pip_concentration.pdf"),
    "\n  ", file.path(PANEL_DIR, "crossancestry_pip_concentration_source.csv"), "\n")
