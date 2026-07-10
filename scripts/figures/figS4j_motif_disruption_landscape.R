#!/usr/bin/env Rscript
# ==============================================================================
# figS4j — Sequence-level TF-motif disruption at PRIORITIZED MASLD target loci.
#
# The sequence-mechanism companion to Fig 4d (which shows open-peak breadth of
# prioritized targets). SAME gene set as every other Fig 4 panel: the canonical
# prioritized universe (fig4a/4d; genetic 3,038 ∪ transcriptomic 8,088 = 9,882),
# protein-coding + the fig4d artifact-symbol filter. Earlier drafts coloured by a
# bespoke "disease master-regulator" TF set (04c: a hepatocyte SCENIC+ regulon TF
# that is ALSO a bulk DEG or COLOC hit) — that set kept 72% of the regulon
# universe and leaned on the TF's OWN mRNA (a poor proxy for TF activity), so it
# is NOT used here; we rely on the paper's canonical prioritization instead.
#
# WHAT IT SHOWS (and does NOT claim): every point is one fine-mapped GWAS
# credible-set variant × TF motif it perturbs (motifbreakR, uniform background).
# `alleleDiff` is the SIGNED PWM-score difference (alt vs ref): >0 = alt/risk
# allele STRENGTHENS the motif, <0 = weakens it. This is a genome-SEQUENCE
# property independent of cohort size (it sidesteps the n=18 multiome
# underpowering) and CORROBORATES the well-powered GWAS/COLOC (Fig 2) + bulk
# RNA-seq (Fig 3). motifbreakR calls are PUTATIVE (not experimentally validated);
# there is NO significance test on the alleleDiff axis (Refpvalue/Altpvalue and
# fimo_concordant are unpopulated in the source table and are not used).
#
#   x    = signed alleleDiff (loss ← 0 → gain)
#   y    = variant fine-map PIP
#   grey = all motif events whose variant does NOT tag a prioritized target
#   ring = the variant tags a PRIORITIZED MASLD target gene (fig4a/4d universe);
#          FILL = that cis-gene's eQTL COLOC PP.H4 (gene-level; same gene as the
#          label, so no trans-TF/cis-gene mismatch), NA (not eQTL-tested) = grey
#   size = motifbreakR effect tier (strong/weak)
#   label= the prioritized cis-GENE (as in fig4a), top hits by priority
#
# Source: GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv
#         + gwas_atac_variant_annotation.csv (variant→cis-gene)
#         + susie_coloc/gene_level_coloc.csv (gene-level COLOC PP.H4)
#         + prioritized_universe_FINAL.txt (fig4a/4d)
# Output: figures/main/fig4_validation/panels/figS4j.pdf   (PDF only)
# Env:    rnaseq  (pure plotting from CSVs — no motifbreakR rerun)
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(ggrastr)
})
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# fig4d's artifact-symbol drop list (non-canonical / artifact-looking; gene-level).
DROP <- c("SPECC1L-ADORA2A","RPS24","EXOC3L4","LRRC75B","SLC50A1","CKM","DYNLRB2","MAU2",
          "GATAD2A","RFXANK","SPPL3","FLAD1","ZNF827","SPATA31H1","SIPA1L2","CWF19L1","FCHO2")

# ── Load motifbreakR results (dedup exact-duplicate rows in the source) ───────
d <- fread(file.path(BASE, "GWAS/finemapping/results/gwas_atac/motif_disruption_scores.csv"))
n_raw <- nrow(d); d <- unique(d)
cat(sprintf("[figS4j] deduplicated source: %d -> %d rows\n", n_raw, nrow(d)))
d <- d[!is.na(alleleDiff) & !is.na(max_pip)]
d[, priority_score := max_pip * abs(alleleDiff)]
d[, effect_tier := factor(ifelse(effect == "strong", "strong", "weak"),
                          levels = c("strong", "weak"))]

# ── Map each variant to its cis-gene, flag prioritized-target membership ──────
va <- fread(file.path(BASE, "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv"))
vg <- va[, .(linked = linked_gene[1], nearest = nearest_gene[1]), by = variant_id]
vg[, cis_gene := fifelse(!is.na(linked) & linked != "", linked, nearest)]
d[, cis_gene := vg$cis_gene[match(SNP_id, vg$variant_id)]]

uni <- trimws(readLines(file.path(BASE, "Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt")))
uni <- uni[uni != ""]
bt  <- fread(file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz"))
pc  <- bt[gene_biotype == "protein_coding", unique(gene_name)]
d[, is_target := !is.na(cis_gene) & cis_gene %in% uni & cis_gene %in% pc & !(cis_gene %in% DROP)]

# Gene-level eQTL COLOC PP.H4 of the cis-gene (same gene as the label).
gl <- fread(file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))
gcol <- intersect(c("gene","symbol","human_symbol"), names(gl))[1]
gl <- gl[, .(cis_coloc = suppressWarnings(max(as.numeric(coloc_best_pp4), na.rm = TRUE))), by = c(gcol)]
gl[is.infinite(cis_coloc), cis_coloc := NA_real_]
d[, cis_coloc := gl$cis_coloc[match(cis_gene, gl[[gcol]])]]

tgt <- d[is_target == TRUE]
cat(sprintf("[figS4j] %d events | %d variants | prioritized-target events = %d (%d variants, %d genes)\n",
            nrow(d), uniqueN(d$SNP_id), nrow(tgt), uniqueN(tgt$SNP_id), uniqueN(tgt$cis_gene)))
cat(sprintf("[figS4j] target cis-genes with eQTL COLOC PP.H4>0.5 = %d | tested = %d | untested(NA) = %d genes\n",
            uniqueN(tgt[cis_coloc > 0.5]$cis_gene), uniqueN(tgt[!is.na(cis_coloc)]$cis_gene),
            uniqueN(tgt[is.na(cis_coloc)]$cis_gene)))

# ── Labels: one per prioritized cis-gene, at its highest-priority event ───────
# Kept deliberately sparse (top N by priority) to avoid crowding the compact panel.
N_LABEL <- 9L
lab <- tgt[order(-priority_score)][!duplicated(cis_gene)][seq_len(min(N_LABEL, uniqueN(tgt$cis_gene)))]

# ── Plot ─────────────────────────────────────────────────────────────────────
xlim <- max(abs(d$alleleDiff)) * 1.02
p <- ggplot() +
  ggrastr::rasterise(
    geom_point(data = d[is_target == FALSE],
               aes(alleleDiff, max_pip, size = effect_tier),
               colour = "#D9D9D9", alpha = 0.45, stroke = 0),
    dpi = 450) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = house_ink) +
  geom_point(data = d[is_target == TRUE],
             aes(alleleDiff, max_pip, fill = cis_coloc, size = effect_tier),
             shape = 21, colour = house_ink, stroke = 0.2) +
  ggrepel::geom_text_repel(
    data = lab, aes(alleleDiff, max_pip, label = cis_gene),
    size = 5.5 / .pt, fontface = "italic", colour = house_ink,
    max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2,
    segment.colour = "grey65", box.padding = 0.9, point.padding = 0.4,
    force = 12, force_pull = 0.4, max.iter = 20000, seed = 7) +
  scale_fill_gradient(low = "#EFEDF5", high = "#54278F", na.value = "#DEDEDE",
                      limits = c(0, 1), breaks = c(0, 0.5, 1), name = "cis-gene eQTL\nCOLOC PP.H4",
                      guide = guide_colourbar(barwidth = unit(0.15, "cm"), barheight = unit(1.0, "cm"), order = 1)) +
  scale_size_manual(values = c(strong = 0.85, weak = 0.32), name = "motifbreakR effect",
                    guide = guide_legend(override.aes = list(shape = 21, fill = "grey50"), order = 2)) +
  scale_x_continuous(limits = c(-xlim, xlim), breaks = c(-2, -1, 0, 1, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.02, 1.04), breaks = c(0, 0.5, 1.0), expand = c(0, 0)) +
  annotate("text", x = -xlim * 0.98, y = 1.03, hjust = 0, vjust = 1,
           label = "motif LOSS (alt weakens)", size = 5.5 / .pt, colour = "grey45") +
  annotate("text", x = xlim * 0.98, y = 1.03, hjust = 1, vjust = 1,
           label = "motif GAIN (alt strengthens)", size = 5.5 / .pt, colour = "grey45") +
  labs(x = "motifbreakR allelic PWM-score difference (alt − ref)",
       y = "variant fine-map PIP") +
  theme_masld() +
  theme(legend.position = "right", legend.key.size = unit(0.28, "lines"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        legend.spacing.y = unit(1, "pt"), legend.margin = margin(0, 0, 0, 2),
        panel.grid.minor = element_blank())

out <- file.path(FIG4_DIR, "panels", "figS4j.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, p, width = 3.7, height = 2.8, device = grDevices::cairo_pdf)
cat("[figS4j] saved:", out, "\n")

# Provenance: the labelled prioritized-target hits (with disrupted TF + coloc).
fwrite(tgt[order(-priority_score)][!duplicated(paste(cis_gene, SNP_id)),
           .(cis_gene, tf_name, SNP_id, alleleDiff, max_pip, priority_score, effect, cis_coloc)],
       file.path(FIG4_DIR, "panels", "figS4j_labeled_hits.csv"))

message(sprintf(paste0(
  "CAPTION (Fig S4j): Sequence-level TF-motif disruption at prioritized MASLD target loci. Each point is ",
  "one fine-mapped GWAS credible-set variant × TF motif event from motifbreakR (%d events, %d variants, ",
  "uniform background). x = signed allelic PWM-score difference (alleleDiff; <0 = risk allele weakens the ",
  "motif, >0 = strengthens); y = variant fine-mapping posterior (max PIP). Grey = motif events whose variant ",
  "does not tag a prioritized target. Ringed points = the %d events (%d variants, %d genes) whose fine-mapped ",
  "variant maps to a gene in the canonical prioritized universe shared across all Fig 4 panels (fig4a/4d; ",
  "genetic 3,038 ∪ transcriptomic 8,088, protein-coding, artifact-filtered), FILL-coloured by that cis-gene's ",
  "gene-level eQTL colocalization (COLOC PP.H4; light-grey = not eQTL-tested). Labels name the prioritized ",
  "cis-gene (the disrupted-motif TF is retained in the source table, not shown, as it is the mechanism not the ",
  "selection). Point size = motifbreakR effect tier (strong/weak). alleleDiff is a genome-sequence property ",
  "independent of cohort size; it is the sequence-mechanism companion to the open-peak breadth of Fig 4d. ",
  "motifbreakR disruptions are PUTATIVE (not experimentally validated); no significance test is applied to the ",
  "alleleDiff axis."),
  nrow(d), uniqueN(d$SNP_id), nrow(tgt), uniqueN(tgt$SNP_id), uniqueN(tgt$cis_gene)))
