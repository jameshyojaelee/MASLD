#!/usr/bin/env Rscript
# figS09_locus_zoom.R
# Locus-zoom visualization for shared EUR+EAS fine-mapping loci.
#
# Generates two complementary plots:
#   1. 4-panel scatter: EUR GWAS / EAS GWAS / PIP EUR / PIP EAS
#      (matches coworker's PDF style — plain gray points, lead SNP highlighted)
#   2. LD-zoom: LD-colored EUR scatter + rotated triangular LD heatmap +
#      finemapped-variant track + eQTL track + gene-annotation track
#      (matches coworker's PNG style)
#
# Usage:
#   Rscript figS09_locus_zoom.R <locus_id> [trait_pair]
#   Rscript figS09_locus_zoom.R locus_ALT_chr19_41353107 ALT
#
# Defaults to locus_ALT_chr19_41353107 (the ALDH2/CLC locus on chr19).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(Matrix)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggrepel)
})

# ── Args ──────────────────────────────────────────────────────────────────────
# Usage: Rscript figS09_locus_zoom.R <locus_id> [trait_pair] [eqtl_gene] [out_ld_zoom_override]
#   arg 4 (optional): if given, the LD-zoom plot writes directly to that absolute
#   path instead of OUT_DIR/{LOCUS_ID}_ld_zoom.pdf — used by the fig3b RORA wrapper.
args            <- commandArgs(trailingOnly = TRUE)
LOCUS_ID        <- if (length(args) >= 1) args[1] else "locus_ALT_chr19_41353107"
TRAIT_PAIR      <- if (length(args) >= 2) args[2] else sub("locus_([^_]+)_.*", "\\1", LOCUS_ID)
FORCE_EQTL_GENE <- if (length(args) >= 3) args[3] else NA_character_  # e.g. "RORA"
OUT_LD_ZOOM_OVERRIDE <- if (length(args) >= 4) args[4] else NA_character_

cat("Locus:", LOCUS_ID, "| Trait:", TRAIT_PAIR, "\n")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
PROJ  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
source(file.path(PROJ, "scripts/figures/publication_theme.R"))

EUR_SS      <- file.path(BASE, "data/sumstats", paste0("UKBB_",  TRAIT_PAIR, "_reformatted_hg19.tsv"))
EAS_SS      <- file.path(BASE, "data/sumstats", paste0("BBJ_",   TRAIT_PAIR, "_reformatted_hg19.tsv"))
MES_VARS    <- file.path(BASE, "results/mesusie/mesusie_variant_summary.csv")
SHARED_LOCI <- file.path(BASE, "results/susiex/shared_loci.csv")
EQTL_DIR    <- file.path(PROJ, "data/broadaway_eqtl")
# LD panel dispatch (default 1kg; sghatan UKBB EUR is off-limits for new computation).
source(file.path(BASE, "src/finemapping_functions.R"))
LD_BASE     <- get_ld_base_dir("EUR")
LD_BLOCKS   <- file.path(LD_BASE, "approx_LD_blocks.txt")
ATLAS_FILE  <- file.path(PROJ, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

OUT_DIR <- file.path(BASE, "figures/locus_zoom")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Load locus metadata ───────────────────────────────────────────────────────
cat("Loading locus metadata...\n")
loci   <- fread(SHARED_LOCI)
locus  <- loci[locus_id == LOCUS_ID]
if (nrow(locus) == 0) stop("Locus not found in shared_loci.csv: ", LOCUS_ID)

CHR          <- locus$chr
WIN_START    <- locus$window_start
WIN_END      <- locus$window_end
EUR_LEAD_POS <- locus$eur_lead_pos
EAS_LEAD_POS <- locus$eas_lead_pos
EUR_LEAD_SNP <- locus$eur_lead_snp   # format: "chr:pos:a1:a2"
EAS_LEAD_SNP <- locus$eas_lead_snp

cat(sprintf("  chr%s:%s-%s | EUR lead: %s | EAS lead: %s\n",
            CHR, WIN_START, WIN_END, EUR_LEAD_POS, EAS_LEAD_POS))

# ── Load MESuSiE PIPs ─────────────────────────────────────────────────────────
cat("Loading MESuSiE PIPs...\n")
pips <- fread(MES_VARS)[locus_id == LOCUS_ID]
cat("  Variants:", nrow(pips), "| In-CS:", sum(pips$in_cs), "\n")

# ── Load EUR sumstats (region-filtered) ───────────────────────────────────────
cat("Loading EUR sumstats (~326MB, filtering to locus window)...\n")
eur_ss <- fread(EUR_SS,
                select = c("chromosome","position","allele1","allele2","beta","se","pval"))
eur_locus <- eur_ss[chromosome == CHR & position >= WIN_START & position <= WIN_END]
eur_locus[, pval    := as.numeric(pval)]
eur_locus[, mlog10p := -log10(pmax(pval, 1e-300))]
eur_locus[, snp_id  := paste(chromosome, position, allele1, allele2, sep = ":")]
cat("  EUR locus variants:", nrow(eur_locus), "\n")
rm(eur_ss); gc()

# ── Load EAS sumstats ──────────────────────────────────────────────────────────
cat("Loading EAS sumstats...\n")
eas_ss <- fread(EAS_SS,
                select = c("chromosome","position","allele1","allele2","beta","se","pval"))
eas_locus <- eas_ss[chromosome == CHR & position >= WIN_START & position <= WIN_END]
eas_locus[, pval    := as.numeric(pval)]
eas_locus[, mlog10p := -log10(pmax(pval, 1e-300))]
cat("  EAS locus variants:", nrow(eas_locus), "\n")
rm(eas_ss); gc()

# ── Merge PIPs into EUR table ─────────────────────────────────────────────────
pips[, snp_id  := paste(chr, pos, a1, a2, sep = ":")]
pips[, snp_id2 := paste(chr, pos, a2, a1, sep = ":")]   # flipped alleles
eur_locus <- merge(
  eur_locus,
  pips[, .(snp_id, pip_eur, pip_eas, cs_category, in_cs)],
  by = "snp_id", all.x = TRUE
)
# Also try flipped match for unmatched
unmatched <- eur_locus[is.na(pip_eur), snp_id]
flip_match <- pips[snp_id2 %in% unmatched, .(snp_id2, pip_eur, pip_eas, cs_category, in_cs)]
if (nrow(flip_match) > 0) {
  eur_locus <- merge(eur_locus, flip_match,
                     by.x = "snp_id", by.y = "snp_id2", all.x = TRUE,
                     suffixes = c("", ".flip"))
  eur_locus[is.na(pip_eur) & !is.na(pip_eur.flip), `:=`(
    pip_eur = pip_eur.flip, pip_eas = pip_eas.flip,
    cs_category = cs_category.flip, in_cs = in_cs.flip
  )]
  eur_locus[, c("pip_eur.flip","pip_eas.flip","cs_category.flip","in_cs.flip") := NULL]
}
cat("  Variants with PIP:", sum(!is.na(eur_locus$pip_eur)), "\n")

# ── Colors ────────────────────────────────────────────────────────────────────
COL_EUR  <- "black"     # Coworker's original black for EUR
COL_EAS  <- "#D55E00"   # Coworker's burnt orange for EAS
GWS      <- -log10(5e-8)

theme_locus <- theme_masld(base_size = 7) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x   = element_blank(),
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.line.y   = element_line(linewidth = 0.5, color = "black"),
    legend.position = "none",
    plot.title    = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5, face = "bold", size = 8, margin = margin(r = 5)),
    strip.background  = element_blank(),
    strip.placement   = "outside",
    plot.margin       = margin(t = 2, r = 15, b = 2, l = 5)
  )
theme_locus_x <- theme_locus +
  theme(
    axis.line.x   = element_line(linewidth = 0.5, color = "black"),
    axis.title.x  = element_text(margin = margin(t = 2)),
    axis.text.x   = element_text(),
    axis.ticks.x  = element_line(linewidth = 0.5, color = "black"),
    plot.margin   = margin(t = 2, r = 15, b = 2, l = 5)
  )

# ════════════════════════════════════════════════════════════════════════════════
# PLOT 1 — 4-panel: EUR / EAS GWAS + PIP EUR / PIP EAS
# ════════════════════════════════════════════════════════════════════════════════
cat("\n── Building Plot 1: 4-panel locus zoom ──\n")

p1_eur <- ggplot(eur_locus, aes(x = position, y = mlog10p)) +
  geom_hline(yintercept = GWS, linetype = "longdash", color = "black", linewidth = 0.4) +
  rasterize_layer(geom_point(color = "grey85", size = 1.2, alpha = 0.4), dpi = 600) +
  geom_point(data = eur_locus[position == EUR_LEAD_POS],
             fill = COL_EUR, color = "black", size = 3.0, shape = 21, stroke = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_x_continuous(limits = c(WIN_START, WIN_END), expand = c(0, 0)) +
  labs(y = expression(-log[10](italic(p))), title = NULL) +
  facet_grid("EUR\nGWAS" ~ ., switch = "y") +
  theme_locus

p1_eas <- ggplot(eas_locus, aes(x = position, y = mlog10p)) +
  geom_hline(yintercept = GWS, linetype = "longdash", color = "black", linewidth = 0.4) +
  rasterize_layer(geom_point(color = "grey85", size = 1.2, alpha = 0.4), dpi = 600) +
  geom_point(data = eas_locus[position == EAS_LEAD_POS],
             fill = COL_EAS, color = "black", size = 3.0, shape = 21, stroke = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_x_continuous(limits = c(WIN_START, WIN_END), expand = c(0, 0)) +
  labs(y = expression(-log[10](italic(p))), title = NULL) +
  facet_grid("EAS\nGWAS" ~ ., switch = "y") +
  theme_locus

p1_pip_eur <- ggplot(pips, aes(x = pos, y = pip_eur)) +
  rasterize_layer(geom_point(color = "grey85", size = 1.2, alpha = 0.4), dpi = 600) +
  geom_point(data = pips[pos == EUR_LEAD_POS],
             fill = COL_EUR, color = "black", size = 3.0, shape = 21, stroke = 0.5) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(WIN_START, WIN_END), expand = c(0, 0)) +
  labs(y = "PIP", title = NULL) +
  facet_grid("PIP\nEUR" ~ ., switch = "y") +
  theme_locus

p1_pip_eas <- ggplot(pips, aes(x = pos, y = pip_eas)) +
  rasterize_layer(geom_point(color = "grey85", size = 1.2, alpha = 0.4), dpi = 600) +
  geom_point(data = pips[pos == EAS_LEAD_POS],
             fill = COL_EAS, color = "black", size = 3.0, shape = 21, stroke = 0.5) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(WIN_START, WIN_END), expand = c(0, 0), labels = function(x) paste0(format(x, big.mark = ",", trim = TRUE), " bp")) +
  labs(x = "Genomic position (hg19)", y = "PIP", title = NULL) +
  facet_grid("PIP\nEAS" ~ ., switch = "y") +
  theme_locus_x

plot1 <- (p1_eur / p1_eas / p1_pip_eur / p1_pip_eas) +
  plot_layout(heights = c(3, 3, 2, 2))

out1 <- file.path(OUT_DIR, paste0(LOCUS_ID, "_4panel.pdf"))
save_fig(plot1, out1, width = fig_half_width, height = 6.5, dpi = 600)
cat("  Saved:", out1, "\n")

# ════════════════════════════════════════════════════════════════════════════════
# PLOT 2 — LD-colored zoom + triangular LD heatmap + tracks
# ════════════════════════════════════════════════════════════════════════════════
cat("\n── Building Plot 2: LD-zoom ──\n")

# ── Find LD block(s) for locus ────────────────────────────────────────────────
blocks_tbl <- fread(LD_BLOCKS)
matching   <- blocks_tbl[chr == CHR & start < WIN_END & stop > WIN_START]
if (nrow(matching) == 0) {
  warning("No LD blocks found — skipping Plot 2"); quit(save = "no")
}
cat("  LD blocks:", nrow(matching), "\n")

# ── Load bim(s) and identify locus row indices ────────────────────────────────
bim_list <- lapply(seq_len(nrow(matching)), function(k) {
  blk <- file.path(LD_BASE, paste0("chr", CHR),
                   paste0(matching$start[k], ".", matching$stop[k]),
                   paste0(matching$start[k], ".", matching$stop[k]))
  bim_file <- paste0(blk, ".bim")
  if (!file.exists(bim_file)) { warning("bim not found: ", bim_file); return(NULL) }
  bim <- fread(bim_file, col.names = c("chr","rsid","dk","pos","alt","ref"))
  bim[, block_path := blk]
  bim[, row_in_block := .I]
  bim
})
bim_all <- rbindlist(bim_list[!sapply(bim_list, is.null)])
bim_locus <- bim_all[pos >= WIN_START & pos <= WIN_END]
cat("  Variants in LD panel within window:", nrow(bim_locus), "\n")

if (nrow(bim_locus) < 10) {
  warning("Too few LD variants in locus window — skipping Plot 2"); quit(save = "no")
}

# ── Load LD matrix (per block), extract locus rows/cols ──────────────────────
cat("  Loading LD matrix (large file — may take ~30s)...\n")
ld_chunks <- lapply(unique(bim_locus$block_path), function(blk_path) {
  ld_file  <- paste0(blk_path, ".ld")
  bim_file <- paste0(blk_path, ".bim")
  if (!file.exists(ld_file)) { warning("LD file not found: ", ld_file); return(NULL) }

  bim_blk  <- fread(bim_file, col.names = c("chr","rsid","dk","pos","alt","ref"))
  idx      <- which(bim_blk$pos >= WIN_START & bim_blk$pos <= WIN_END)
  if (length(idx) == 0) return(NULL)

  # Load full LD block and subset immediately
  ld_full  <- as.matrix(fread(ld_file))
  ld_sub   <- ld_full[idx, idx, drop = FALSE]
  colnames(ld_sub) <- rownames(ld_sub) <-
    paste(bim_blk$chr[idx], bim_blk$pos[idx], bim_blk$alt[idx], bim_blk$ref[idx], sep = ":")
  rm(ld_full); gc()
  list(ld = ld_sub, bim = bim_blk[idx])
})
ld_chunks <- ld_chunks[!sapply(ld_chunks, is.null)]

# Combine blocks (block-diagonal if multiple)
if (length(ld_chunks) == 1) {
  ld_mat  <- ld_chunks[[1]]$ld
  bim_sub <- ld_chunks[[1]]$bim
} else {
  bim_sub <- rbindlist(lapply(ld_chunks, `[[`, "bim"))
  mats    <- lapply(ld_chunks, `[[`, "ld")
  ld_mat  <- as.matrix(bdiag(mats))
  snp_ids <- unlist(lapply(mats, rownames))
  colnames(ld_mat) <- rownames(ld_mat) <- snp_ids
}
cat("  LD submatrix:", nrow(ld_mat), "×", ncol(ld_mat), "\n")

# ── r² to EUR lead SNP ────────────────────────────────────────────────────────
# Match by position (allele ordering differs between LD panel and sumstats)
ld_pos_vec <- as.integer(sub("^[^:]+:([^:]+):.*", "\\1", rownames(ld_mat)))
lead_row   <- which.min(abs(ld_pos_vec - EUR_LEAD_POS))
cat("  EUR lead matched at LD pos", ld_pos_vec[lead_row],
    "(target:", EUR_LEAD_POS, ")\n")

r_to_lead  <- ld_mat[lead_row[1], ]   # correlation (r)
r2_to_lead <- as.numeric(r_to_lead)^2 # r²
r2_by_pos  <- data.table(pos = ld_pos_vec, r2 = r2_to_lead)
# Keep max r² per position (handles multi-allelic sites)
r2_by_pos  <- r2_by_pos[, .(r2 = max(r2, na.rm = TRUE)), by = pos]
r2_by_pos[r2 > 1, r2 := 1]

# Merge r² into EUR locus by position (r2_by_pos uses "pos"; eur_locus uses "position")
setnames(r2_by_pos, "pos", "position")
eur_locus <- merge(eur_locus, r2_by_pos, by = "position", all.x = TRUE)
# Discretize r2 for scatter
eur_locus[, r2_class := cut(r2, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, 1.0),
                            labels = c("0.0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1.0"))]
# Ensure NAs are handled (they stay NA, which ggplot maps correctly to missing)

# ── LD color scales ───────────────────────────────────────────────────────────
# Scatter: near-white → cyan → violet → magenta.
# Low end uses fibrosis_stage_colors$F0 ("#E3F2FD") from publication_theme.R —
# near-white with the faintest blue tint, distinct from white background.
# High end stays on the theme's magenta axis (#E91E63 / #C2185B).
r2_discrete_pal <- c(
  "0.0 - 0.2" = "#E3F2FD",   # near-white light blue (fibrosis F0, publication theme)
  "0.2 - 0.4" = "#00BCD4",   # cyan
  "0.4 - 0.6" = "#7B1FA2",   # violet (masld_colors$human_enriched)
  "0.6 - 0.8" = "#E91E63",   # bright magenta (masld_colors$twas)
  "0.8 - 1.0" = "#C2185B",   # deep magenta (masld_colors$up)
  "Lead SNP"  = "#212121"    # near-black diamond
)

# LD triangle heatmap: YlGnBu (ColorBrewer) — near-white cream → light green →
# cyan → deep navy. Classic LD heatmap scale; matches "dark blue/cyan to
# yellowish/light green almost white" description exactly.
r2_pal_heatmap <- c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB",
                    "#41B6C4", "#1D91C0", "#225EA8", "#081D58")

# ── LD-colored EUR scatter ───────────────────────────────────────────────────
# Sort so high-r² points plot on top
eur_locus <- eur_locus[order(r2)]
lead_pt   <- eur_locus[position == EUR_LEAD_POS]
# Add "Lead SNP" as fill class so it appears in legend
lead_pt[, r2_class := factor("Lead SNP",
  levels = c(levels(eur_locus$r2_class), "Lead SNP"))]
eur_locus[, r2_class := factor(r2_class,
  levels = c(levels(r2_class), "Lead SNP"))]

p2_eur_ld <- ggplot(eur_locus, aes(x = position, y = mlog10p)) +
  geom_hline(yintercept = GWS, linetype = "longdash", color = "black", linewidth = 0.4) +
  rasterize_layer(geom_point(aes(fill = r2_class), color = "transparent", shape = 21, size = 1.8, alpha = 0.85), dpi = 600) +
  geom_point(data = lead_pt, aes(fill = r2_class), color = "black", size = 1.5, shape = 23, stroke = 0.5) +
  scale_fill_manual(
    values = r2_discrete_pal,
    na.value = "grey85",
    name = expression(r^2),
    breaks = c("0.0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1.0", "Lead SNP"),
    guide = guide_legend(
      override.aes = list(
        shape  = c(21, 21, 21, 21, 21, 23),
        color  = c(rep(NA, 5), "black"),
        size   = c(rep(2.5, 5), 2.5)
      )
    )
  ) +
  scale_x_continuous(limits = c(WIN_START, WIN_END), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(y = expression(-log[10](italic(p))), title = NULL) +
  facet_grid("Summary\nstatistics" ~ ., switch = "y") +
  theme_locus +
  theme(legend.position = "right", legend.title = element_text(size = 8),
        legend.key.height = unit(0.3, "cm"), legend.key = element_rect(color="white"),
        plot.margin = margin(b = 0, unit = "pt"))

# ── Rotated triangular LD heatmap ─────────────────────────────────────────────
# Standard LocusZoom approach:
#   For each pair (i, j) with i <= j (upper triangle):
#     x_center = (pos_i + pos_j) / 2
#     y_depth  = -(pos_j - pos_i) / 2   (below x-axis)
#   Plot as geom_tile with equal width/height = mean variant spacing

  # Thin LD variants uniformly from the LD panel (NOT sumstat-matched).
  # Use ~400 variants: fewer pairs (~80K vs 1.3M), larger effective spacing,
  # and tiles that fill the triangle visually at typical figure heights.
  N_target    <- 400
  thin_step   <- max(1L, floor(length(ld_pos_vec) / N_target))
  keep_ld_idx <- seq(1L, length(ld_pos_vec), by = thin_step)
  ld_mat_thin   <- ld_mat[keep_ld_idx, keep_ld_idx, drop = FALSE]
  snp_pos_thin  <- ld_pos_vec[keep_ld_idx]
  snp_pos_ord   <- sort(unique(snp_pos_thin))
  n_snps        <- length(snp_pos_ord)
  cat("  LD heatmap variants after thinning:", n_snps, "\n")

if (n_snps >= 2) {
  cat("  Building triangular LD heatmap...\n")

  # Map snp names to sorted indices
  pos2idx      <- setNames(seq_along(snp_pos_ord), as.character(snp_pos_ord))

  # Build upper-triangle pair table
  pairs <- data.table(
    i = rep(seq_len(n_snps), each = n_snps),
    j = rep(seq_len(n_snps), n_snps)
  )[i <= j]
  pairs[, p_i := snp_pos_ord[i]]
  pairs[, p_j := snp_pos_ord[j]]

  # Look up r² for each pair
  # Build pos→row index from thinned LD matrix
  pos_to_ldrow <- setNames(seq_along(snp_pos_thin), as.character(snp_pos_thin))
  pairs[, row_i := pos_to_ldrow[as.character(p_i)]]
  pairs[, row_j := pos_to_ldrow[as.character(p_j)]]
  valid <- !is.na(pairs$row_i) & !is.na(pairs$row_j)
  pairs <- pairs[valid]

  if (nrow(pairs) > 0) {
    r2_vals <- ld_mat_thin[cbind(pairs$row_i, pairs$row_j)]^2
    pairs[, r2_pair := r2_vals]
    pairs[r2_pair > 1, r2_pair := 1]
    # Pairs with r²<0.01 are set to NA and rendered as light gray (#E8E8E8),
    # making the triangular background visible against the white panel.

    # Rotated coordinates
    pairs[, x_center := (p_i + p_j) / 2]
    pairs[, y_depth  := -(p_j - p_i) / 2]
    # Larger tiles (window/30) fill the visible region more densely now that
    # the triangle is heavily truncated.
    avg_sp        <- mean(diff(snp_pos_ord))
    min_tile_size <- (WIN_END - WIN_START) / 30
    tile_size     <- max(avg_sp * max(thin_step, 8), min_tile_size)

    # Show only top 20% of max depth — keeps the informative near-diagonal
    # LD blocks and drops the distant low-r² tail entirely.
    tri_y_min <- -(WIN_END - WIN_START) / 2 * 0.20

    p2_ld_tri <- ggplot(pairs, aes(x = x_center, y = y_depth, fill = r2_pair)) +
      rasterize_layer(geom_tile(width = tile_size, height = tile_size), dpi = 600) +
      scale_fill_gradientn(
        colors = r2_pal_heatmap, limits = c(0, 1), na.value = "#E8E8E8",
        name = expression(r^2)
      ) +
      scale_x_continuous(limits = c(WIN_START, WIN_END), expand = c(0, 0)) +
      scale_y_continuous(limits = c(tri_y_min, 0), expand = c(0, 0)) +
      labs(y = "dummy", title = NULL) +
      facet_grid("LD matrix" ~ ., switch = "y") +
      theme_locus +
      theme(
        legend.position = "none",
        axis.title.y = element_text(color = "transparent"),
        axis.text.y = element_text(color = "transparent"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = margin(t = 0, unit = "pt")
      )
  } else {
    p2_ld_tri <- ggplot() + annotate("text", x=0.5, y=0.5, label="No LD pairs") + theme_void()
  }
} else {
  cat("  No LD heatmap — fewer than 2 thinned variants\n")
  p2_ld_tri <- ggplot() + theme_void()
}

# ── Finemapped variants track ─────────────────────────────────────────────────
# Color by CS index (CS1, CS2, ...) — more informative than ancestry category
in_cs_vars <- pips[in_cs == TRUE]

# Assign a per-locus CS label: combine ancestry + index for distinct credible sets
in_cs_vars[, cs_label := paste0(cs_category, " CS", cs_index)]

# CS colors: EUR=reds, EAS=blues, Shared=purples, then cycle
# Use ancestry-aware colors so EUR/EAS/Shared are immediately distinguishable
cs_labels_ordered <- unique(in_cs_vars[order(cs_category, cs_index), cs_label])
cs_pal_map <- function(lab) {
  if      (grepl("^EUR",    lab)) c("#D73027","#FC8D59","#FDAE61","#FEE090")[1]
  else if (grepl("^EAS",    lab)) c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8")[1]
  else if (grepl("^Shared", lab)) c("#762A83","#C2A5CF","#A6D96A","#1A9850")[1]
  else                             "#AAAAAA"
}
# Assign slightly different shades for multiple CS per ancestry
assign_cs_color <- function(labs) {
  eur_pal    <- c("#D73027","#FC8D59","#FDAE61","#FEE090")
  eas_pal    <- c("#4575B4","#74ADD1","#ABD9E9")
  shared_pal <- c("#762A83","#C2A5CF","#1A9850")
  eur_i <- eas_i <- shared_i <- 1L
  setNames(sapply(labs, function(lab) {
    if      (grepl("^EUR",    lab)) { col <- eur_pal[min(eur_i,    length(eur_pal))];    eur_i    <<- eur_i    + 1L; col }
    else if (grepl("^EAS",    lab)) { col <- eas_pal[min(eas_i,    length(eas_pal))];    eas_i    <<- eas_i    + 1L; col }
    else if (grepl("^Shared", lab)) { col <- shared_pal[min(shared_i,length(shared_pal))];shared_i <<- shared_i + 1L; col }
    else "#AAAAAA"
  }), labs)
}
cs_colors <- assign_cs_color(cs_labels_ordered)

p2_cs <- ggplot() +
  # Threshold line at PIP = 0.9
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  # Lollipop stems for all variants
  geom_segment(data = pips, aes(x = pos, xend = pos, y = 0, yend = pip * 0.9),
               color = "grey82", linewidth = 0.35) +
  # All variants (milky grey background)
  rasterize_layer(geom_point(data = pips[in_cs == FALSE],
             aes(x = pos, y = pip), fill = "grey85", color = "transparent", shape = 21, size = 1.8, alpha = 0.5), dpi = 600) +
  # CS variants colored by CS index
  geom_point(data = in_cs_vars,
             aes(x = pos, y = pip, color = cs_label),
             size = 1.8) +
  scale_color_manual(values = cs_colors, name = "Credible set") +
  scale_x_continuous(limits = c(WIN_START, WIN_END), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.15), expand = c(0, 0)) +
  labs(y = "PIP", title = NULL) +
  facet_grid("Finemapped\nvariants" ~ ., switch = "y") +
  theme_locus_x +
  theme(axis.line.x = element_line(linewidth=0.5, color="black"),
        legend.position = "right", legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"))

# ── eQTL track ────────────────────────────────────────────────────────────────
eqtl_file <- file.path(EQTL_DIR, paste0("chr", CHR, "_marginal_summary_results.tsv"))
p2_eqtl   <- NULL

if (file.exists(eqtl_file)) {
  eqtl_chr <- fread(eqtl_file,
                    select = c("Entrez","Variant","CHR","POS","Beta","SE","PVAL","GeneSymbol","ENSG"))
  eqtl_locus <- eqtl_chr[CHR == CHR & POS >= WIN_START & POS <= WIN_END]

  if (nrow(eqtl_locus) > 0) {
    # Gene selection:
    # (a) If FORCE_EQTL_GENE is provided and has data in this window, use it.
    # (b) Otherwise, pick the gene with the most significant eQTL within ±75kb.
    if (!is.na(FORCE_EQTL_GENE) && FORCE_EQTL_GENE %in% eqtl_locus$GeneSymbol) {
      best_gene  <- FORCE_EQTL_GENE
      cat("  eQTL track: using FORCE_EQTL_GENE =", best_gene, "\n")
    } else {
      if (!is.na(FORCE_EQTL_GENE))
        warning("FORCE_EQTL_GENE '", FORCE_EQTL_GENE,
                "' not found in locus eQTL data — falling back to auto-selection.")
      near_lead  <- eqtl_locus[abs(POS - EUR_LEAD_POS) <= 75000]
      if (nrow(near_lead) == 0) near_lead <- eqtl_locus
      best_gene  <- near_lead[which.min(PVAL), GeneSymbol]
    }

    eqtl_top <- eqtl_locus[GeneSymbol == best_gene]
    eqtl_top[, mlog10p_eqtl := -log10(pmax(PVAL, 1e-300))]

    # Label as "(nearby)" when the lead SNP is NOT close to any eQTL variant
    # for best_gene — signals the eQTL gene differs from the GWAS gene of interest.
    lead_in_gene <- nrow(eqtl_locus[GeneSymbol == best_gene &
                                    abs(POS - EUR_LEAD_POS) <= 10000]) > 0
    strip_text <- if (lead_in_gene || !is.na(FORCE_EQTL_GENE)) {
      sprintf("eQTL\n%s", best_gene)   # no (nearby) when gene was forced
    } else {
      sprintf("eQTL\n%s\n(nearby)", best_gene)
    }
    eqtl_top[, eqtl_strip := strip_text]

    p2_eqtl <- ggplot(eqtl_top, aes(x = POS, y = mlog10p_eqtl)) +
      rasterize_layer(geom_point(color = "#E07B39", size = 1.0, alpha = 0.80), dpi = 600) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      labs(y = expression(-log[10](italic(p))), title = NULL) +
      facet_grid(eqtl_strip ~ ., switch = "y") +
      theme_locus
    cat("  eQTL track: gene", best_gene, "(nearby:", !lead_in_gene, ")|",
        nrow(eqtl_top), "variants\n")
  }
}

# ── Gene annotation track (TxDb hg19 — exon-level) ───────────────────────────
cat("  Building gene track from TxDb.Hsapiens.UCSC.hg19.knownGene...\n")
txdb     <- TxDb.Hsapiens.UCSC.hg19.knownGene
locus_gr <- GRanges(seqnames = paste0("chr", CHR),
                    ranges   = IRanges(WIN_START, WIN_END))

# Gene-level metadata (strand, DEG status)
genes_gr <- suppressMessages(genes(txdb, filter = list(tx_chrom = paste0("chr", CHR))))
hits     <- subsetByOverlaps(genes_gr, locus_gr)
genes_df <- data.frame(
  entrez = hits$gene_id,
  g_start = start(hits),
  g_end   = end(hits),
  strand  = as.character(strand(hits)),
  stringsAsFactors = FALSE
)
genes_df$symbol <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = genes_df$entrez,
         column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
)
genes_df <- genes_df[!is.na(genes_df$symbol), ]

# DEG status from Atlas
if (file.exists(ATLAS_FILE)) {
  atlas <- fread(ATLAS_FILE, select = c("human_symbol", "bulk_logFC", "bulk_padj"))
  stopifnot(all(c("bulk_padj", "bulk_logFC") %in% names(atlas)))
  genes_df <- merge(genes_df, atlas, by.x = "symbol", by.y = "human_symbol", all.x = TRUE)
} else {
  genes_df$bulk_logFC <- NA_real_
  genes_df$bulk_padj  <- NA_real_
}
genes_df$deg_status <- "Not Significant"
genes_df$deg_status[!is.na(genes_df$bulk_padj) & genes_df$bulk_padj < 0.05 & genes_df$bulk_logFC >  0.5] <- "Upregulated"
genes_df$deg_status[!is.na(genes_df$bulk_padj) & genes_df$bulk_padj < 0.05 & genes_df$bulk_logFC < -0.5] <- "Downregulated"
genes_df$deg_status <- factor(genes_df$deg_status, levels = c("Upregulated", "Downregulated", "Not Significant"))

# Exon rectangles: for each gene pick the transcript with most exons in window
exons_by_tx <- exonsBy(txdb, by = "tx")
tx_gr_all   <- suppressMessages(transcripts(txdb))
tx_in       <- subsetByOverlaps(tx_gr_all, locus_gr)
tx2gene     <- suppressMessages(
  select(txdb, keys = as.character(tx_in$tx_id),
         columns = c("TXID", "GENEID"), keytype = "TXID")
)
tx2gene$symbol <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = tx2gene$GENEID,
         column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
)
tx2gene <- tx2gene[!is.na(tx2gene$symbol) & tx2gene$symbol %in% genes_df$symbol, ]

exon_df <- do.call(rbind, lapply(unique(tx2gene$symbol), function(sym) {
  txids <- tx2gene$TXID[tx2gene$symbol == sym]
  best  <- NULL; best_n <- -1L
  for (tid in as.character(txids)) {
    if (!tid %in% names(exons_by_tx)) next
    e_in <- subsetByOverlaps(exons_by_tx[[tid]], locus_gr)
    if (length(e_in) > best_n) { best <- e_in; best_n <- length(e_in) }
  }
  if (is.null(best) || best_n == 0) return(NULL)
  data.frame(symbol  = sym,
             ex_start = pmax(start(best), WIN_START),
             ex_end   = pmin(end(best),   WIN_END),
             stringsAsFactors = FALSE)
}))

# Backbone: gene body clipped to window
genes_df$start <- pmax(genes_df$g_start, WIN_START)
genes_df$end   <- pmin(genes_df$g_end,   WIN_END)
genes_df$mid   <- (genes_df$start + genes_df$end) / 2
genes_df        <- genes_df[order(genes_df$start), ]
cat("  Genes in window:", nrow(genes_df), "\n")

# Greedy 2-row stagger
genes_df$row <- 1L
row_end <- c(-Inf, -Inf)
gap     <- (WIN_END - WIN_START) * 0.02
for (i in seq_len(nrow(genes_df))) {
  placed <- FALSE
  for (r in 1:2) {
    if (genes_df$start[i] > row_end[r] + gap) {
      genes_df$row[i] <- r; row_end[r] <- genes_df$end[i]
      placed <- TRUE; break
    }
  }
  if (!placed) genes_df$row[i] <- 1L
}

# Propagate row and DEG status to exon_df
if (!is.null(exon_df) && nrow(exon_df) > 0) {
  exon_df <- merge(exon_df,
                   genes_df[, c("symbol", "row", "deg_status")],
                   by = "symbol", all.x = TRUE)
}

# Strand arrow at the visible end of each gene
genes_df$arrow_x <- ifelse(genes_df$strand == "+", genes_df$end, genes_df$start)

p2_genes <- ggplot(genes_df) +
  # Thin backbone (intron line)
  geom_segment(aes(x = start, xend = end, y = row, yend = row, color = deg_status),
               linewidth = 0.35, lineend = "butt") +
  # Strand arrow
  geom_segment(aes(x = arrow_x,
                   xend = ifelse(strand == "+",
                     pmin(arrow_x + (WIN_END - WIN_START) * 0.012, WIN_END),
                     pmax(arrow_x - (WIN_END - WIN_START) * 0.012, WIN_START)),
                   y = row, yend = row, color = deg_status),
               arrow = arrow(length = unit(2.5, "pt"), type = "closed"),
               linewidth = 0.5) +
  # Exon boxes
  { if (!is.null(exon_df) && nrow(exon_df) > 0)
      geom_rect(data = exon_df,
                aes(xmin = ex_start, xmax = ex_end,
                    ymin = row - 0.18, ymax = row + 0.18,
                    fill = deg_status),
                color = NA)
    else list() } +
  # Gene labels
  geom_text(aes(x = mid, y = row - 0.38, label = symbol),
            size = 2.0, fontface = "italic", hjust = 0.5, color = "grey20") +
  scale_color_manual(values = c("Upregulated" = "#D73027", "Downregulated" = "#4575B4",
                                 "Not Significant" = "grey65"), name = "MASLD DEG") +
  scale_fill_manual(values  = c("Upregulated" = "#D73027", "Downregulated" = "#4575B4",
                                 "Not Significant" = "grey65"), guide = "none") +
  scale_x_continuous(limits = c(WIN_START, WIN_END), expand = c(0, 0),
                     labels = function(x) paste0(round(x / 1e6, 2), " Mb")) +
  scale_y_continuous(limits = c(0.4, 2.8), expand = c(0, 0)) +
  labs(x = "Genomic position (hg19)", y = "dummy", title = NULL) +
  facet_grid("Genes" ~ ., switch = "y") +
  theme_locus_x +
  theme(
    axis.title.y  = element_text(color = "transparent"),
    axis.text.y   = element_text(color = "transparent"),
    axis.ticks.y  = element_blank(),
    axis.line.y   = element_blank(),
    axis.text.x   = element_text(size = 8, color = "grey30"),
    axis.ticks.x  = element_line(color = "grey50"),
    axis.line.x   = element_line(color = "grey50"),
    legend.position    = "right",
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 7),
    legend.key.height  = unit(0.3, "cm"),
    legend.key         = element_rect(fill = "white")
  )

# ── Assemble Plot 2 ───────────────────────────────────────────────────────────
panels <- list(
  p2_eur_ld,   # LD-colored EUR GWAS scatter
  p2_ld_tri,   # Rotated triangular LD heatmap
  p2_cs,       # Finemapped variants
  p2_genes     # Gene annotation
)
heights <- c(3, 2, 1.5, 1.2)

if (!is.null(p2_eqtl)) {
  panels  <- append(panels, list(p2_eqtl), after = 2)
  heights <- c(3, 2, 1.5, 1.5, 1.2)
}

plot2 <- Reduce(`/`, panels) +
  plot_layout(heights = heights)

out2 <- if (!is.na(OUT_LD_ZOOM_OVERRIDE)) {
  dir.create(dirname(OUT_LD_ZOOM_OVERRIDE), recursive = TRUE, showWarnings = FALSE)
  OUT_LD_ZOOM_OVERRIDE
} else {
  file.path(OUT_DIR, paste0(LOCUS_ID, "_ld_zoom.pdf"))
}
save_fig(plot2, out2, width = fig_full_width, height = 7, dpi = 600)
cat("  Saved:", out2, "\n")

cat("\nDone.\n")
cat("  Plot 1:", out1, "\n")
cat("  Plot 2:", out2, "\n")
