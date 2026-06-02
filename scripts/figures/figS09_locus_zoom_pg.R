#!/usr/bin/env Rscript
# figS09_locus_zoom_pg.R — plotgardener rewrite of the LD-zoom locus panel.
#
# Mirrors figS09_locus_zoom.R Plot 2 (LD-zoom) but renders with plotgardener
# tracks for coordinate-perfect alignment and accurate genomic axis ticks.
#
# Tracks (top → bottom, all sharing the same pgParams):
#   1. EUR GWAS Manhattan, LD-colored (r² to EUR lead)
#   2. EAS GWAS Manhattan
#   3. PIP track (EUR + EAS lollipops, custom grid primitives)
#   4. eQTL Manhattan (focal gene; auto-pick or FORCE_EQTL_GENE)
#   5. plotgardener gene track (TxDb hg19, DEG-status highlight)
#   6. plotGenomeLabel ruler (Mb)
#
# Usage:
#   Rscript figS09_locus_zoom_pg.R <locus_id> [trait_pair] [eqtl_gene] [out_override]
# Run inside the `plotgardener` micromamba env (R 4.4.3 + plotgardener 1.12.0).

suppressPackageStartupMessages({
  library(data.table)
  library(plotgardener)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(grid)
  library(Matrix)
})

# ── Args ──────────────────────────────────────────────────────────────────────
args            <- commandArgs(trailingOnly = TRUE)
LOCUS_ID        <- if (length(args) >= 1) args[1] else "locus_GGT_chr15_60883281"
TRAIT_PAIR      <- if (length(args) >= 2) args[2] else sub("locus_([^_]+)_.*", "\\1", LOCUS_ID)
FORCE_EQTL_GENE <- if (length(args) >= 3) args[3] else NA_character_
OUT_OVERRIDE    <- if (length(args) >= 4) args[4] else NA_character_

cat(sprintf("Locus: %s | Trait: %s | eQTL gene: %s\n",
            LOCUS_ID, TRAIT_PAIR,
            if (is.na(FORCE_EQTL_GENE)) "auto" else FORCE_EQTL_GENE))

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
PROJ  <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

EUR_SS      <- file.path(BASE, "data/sumstats", paste0("UKBB_", TRAIT_PAIR, "_reformatted_hg19.tsv"))
EAS_SS      <- file.path(BASE, "data/sumstats", paste0("BBJ_",  TRAIT_PAIR, "_reformatted_hg19.tsv"))
MES_VARS    <- file.path(BASE, "results/mesusie/mesusie_variant_summary.csv")
SHARED_LOCI <- file.path(BASE, "results/susiex/shared_loci.csv")
EQTL_DIR    <- file.path(PROJ, "data/broadaway_eqtl")
ATLAS_FILE  <- file.path(PROJ, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")

# LD panel dispatch (default 1kg)
source(file.path(BASE, "src/finemapping_functions.R"))
LD_BASE   <- get_ld_base_dir("EUR")
LD_BLOCKS <- file.path(LD_BASE, "approx_LD_blocks.txt")
EAS_LD_BASE   <- get_ld_base_dir("EAS")
EAS_LD_BLOCKS <- file.path(EAS_LD_BASE, "approx_LD_blocks.txt")

OUT_DIR <- file.path(BASE, "figures/locus_zoom")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Locus metadata ────────────────────────────────────────────────────────────
loci  <- fread(SHARED_LOCI)
locus <- loci[locus_id == LOCUS_ID]
if (nrow(locus) == 0) stop("Locus not found in shared_loci.csv: ", LOCUS_ID)

CHR          <- locus$chr
WIN_START    <- locus$window_start
WIN_END      <- locus$window_end
EUR_LEAD_POS <- locus$eur_lead_pos
EAS_LEAD_POS <- locus$eas_lead_pos

cat(sprintf("  chr%s:%s-%s | EUR lead: %s | EAS lead: %s\n",
            CHR, WIN_START, WIN_END, EUR_LEAD_POS, EAS_LEAD_POS))

# ── EUR sumstats ──────────────────────────────────────────────────────────────
cat("Loading EUR sumstats...\n")
eur_ss <- fread(EUR_SS,
                select = c("chromosome","position","allele1","allele2","beta","se","pval"))
eur_locus <- eur_ss[chromosome == CHR & position >= WIN_START & position <= WIN_END]
eur_locus[, p := as.numeric(pval)]
eur_locus[, snp := paste(chromosome, position, allele1, allele2, sep = ":")]
rm(eur_ss); gc(verbose = FALSE)
cat("  EUR locus variants:", nrow(eur_locus), "\n")

# ── EAS sumstats ──────────────────────────────────────────────────────────────
cat("Loading EAS sumstats...\n")
eas_ss <- fread(EAS_SS,
                select = c("chromosome","position","allele1","allele2","beta","se","pval"))
eas_locus <- eas_ss[chromosome == CHR & position >= WIN_START & position <= WIN_END]
eas_locus[, p := as.numeric(pval)]
rm(eas_ss); gc(verbose = FALSE)
cat("  EAS locus variants:", nrow(eas_locus), "\n")

# ── Single-ancestry SuSiE finemapping PIPs (replaces prior MeSuSiE track) ─────
# Pulls per-variant SuSiE PIPs from combined_finemapping_1kg.csv (canonical
# post-2026-04-23 for non-UKBB; UKBB EUR studies use PolyFun and don't appear
# in this aggregate). Tries UKBB_<TRAIT> + BBJ_<TRAIT> for liver-enzyme
# panels; otherwise uses TRAIT_PAIR directly as study name.
fm_studies <- if (TRAIT_PAIR %in% c("ALT","AST","GGT")) {
  c(paste0("UKBB_", TRAIT_PAIR), paste0("BBJ_", TRAIT_PAIR))
} else {
  TRAIT_PAIR
}
fm_combined <- fread(file.path(BASE, "results/combined_finemapping_1kg.csv"),
  select = c("study","chromosome","position","susie_pip","susie_cs","susie_converged"))
pips <- fm_combined[study %in% fm_studies & chromosome == CHR &
                    position >= WIN_START & position <= WIN_END &
                    susie_converged == TRUE & !is.na(susie_pip) & susie_pip > 0]
cat(sprintf("  SuSiE finemap PIPs (combined CSV): %d variants across %d studies\n",
            nrow(pips), uniqueN(pips$study)))

# (Inline UKBB EUR SuSiE moved below — needs the LD matrix loaded first.)

# ── LD r² loader (shared between EUR + EAS Manhattans) ───────────────────────
load_ld_for_window <- function(ld_base, ld_blocks_file) {
  blocks_tbl <- fread(ld_blocks_file)
  matching   <- blocks_tbl[chr == CHR & start < WIN_END & stop > WIN_START]
  ld_chunks <- lapply(seq_len(nrow(matching)), function(k) {
    blk <- file.path(ld_base, paste0("chr", CHR),
                     paste0(matching$start[k], ".", matching$stop[k]),
                     paste0(matching$start[k], ".", matching$stop[k]))
    bim_file <- paste0(blk, ".bim")
    ld_file  <- paste0(blk, ".ld")
    if (!file.exists(bim_file) || !file.exists(ld_file)) return(NULL)
    bim <- fread(bim_file, col.names = c("chr","rsid","dk","pos","alt","ref"))
    idx <- which(bim$pos >= WIN_START & bim$pos <= WIN_END)
    if (!length(idx)) return(NULL)
    ld_full <- as.matrix(fread(ld_file))
    ld_sub  <- ld_full[idx, idx, drop = FALSE]
    rownames(ld_sub) <- colnames(ld_sub) <- bim$pos[idx]
    rm(ld_full); gc(verbose = FALSE)
    list(ld = ld_sub, pos = bim$pos[idx])
  })
  ld_chunks <- ld_chunks[!sapply(ld_chunks, is.null)]
  if (!length(ld_chunks)) return(NULL)
  if (length(ld_chunks) == 1) {
    list(mat = ld_chunks[[1]]$ld, pos = ld_chunks[[1]]$pos)
  } else {
    positions <- unlist(lapply(ld_chunks, `[[`, "pos"))
    mat <- as.matrix(Matrix::bdiag(lapply(ld_chunks, `[[`, "ld")))
    rownames(mat) <- colnames(mat) <- positions
    list(mat = mat, pos = positions)
  }
}

r2_to_lead <- function(ld_obj, lead_pos) {
  if (is.null(ld_obj)) return(data.table(position = integer(0), r2 = double(0)))
  lead_row <- which.min(abs(ld_obj$pos - lead_pos))
  v <- as.numeric(ld_obj$mat[lead_row, ])^2
  v[v > 1] <- 1
  dt <- data.table(position = as.integer(ld_obj$pos), r2 = v)
  dt[, .(r2 = max(r2, na.rm = TRUE)), by = position]
}

# EUR LD
cat("Loading EUR LD matrix and computing r² to EUR lead...\n")
eur_ld   <- load_ld_for_window(LD_BASE, LD_BLOCKS)
ld_mat   <- eur_ld$mat
ld_pos   <- eur_ld$pos
eur_r2   <- r2_to_lead(eur_ld, EUR_LEAD_POS)
eur_locus <- merge(eur_locus, eur_r2, by = "position", all.x = TRUE)
eur_locus[is.na(r2), r2 := 0]
cat(sprintf("  LD merged: %d EUR variants matched to EUR LD panel (out of %d)\n",
            sum(!is.na(eur_locus$r2)), nrow(eur_locus)))

# EAS LD (parallel to EUR — uses 1kg_eas panel)
cat("Loading EAS LD matrix and computing r² to EAS lead...\n")
eas_ld <- if (file.exists(EAS_LD_BLOCKS)) {
            load_ld_for_window(EAS_LD_BASE, EAS_LD_BLOCKS)
          } else NULL
if (!is.null(eas_ld)) {
  eas_r2 <- r2_to_lead(eas_ld, EAS_LEAD_POS)
  eas_locus <- merge(eas_locus, eas_r2, by = "position", all.x = TRUE)
  eas_locus[is.na(r2), r2 := 0]
  cat(sprintf("  LD merged: %d EAS variants matched to EAS LD panel (out of %d)\n",
              sum(!is.na(eas_locus$r2)), nrow(eas_locus)))
} else {
  eas_locus[, r2 := 0]
  cat("  No EAS LD blocks file found; EAS panel will render without LD coloring.\n")
}

# ── Inline UKBB EUR SuSiE (only if no UKBB study had chr-level data) ───────
ukbb_present <- any(grepl("^UKBB", unique(pips$study)))
if (!ukbb_present && requireNamespace("susieR", quietly = TRUE)) {
  cat("  Running inline SuSiE on UKBB EUR sumstats (no combined-CSV hit)...\n")
  shared_pos <- intersect(eur_locus$position, as.integer(rownames(ld_mat)))
  if (length(shared_pos) >= 30) {
    eur_sub <- eur_locus[position %in% shared_pos][order(position)]
    pos_chr <- as.character(eur_sub$position)
    ld_sub  <- ld_mat[pos_chr, pos_chr]
    z <- if (all(c("beta","se") %in% names(eur_sub))) {
           eur_sub$beta / eur_sub$se
         } else if ("p" %in% names(eur_sub)) {
           # sign from beta if available, else from positive convention
           sgn <- if ("beta" %in% names(eur_sub)) sign(eur_sub$beta) else 1
           sgn * qnorm(pmin(1 - 1e-300, eur_sub$p / 2), lower.tail = FALSE)
         } else NULL
    if (!is.null(z) && length(z) == nrow(ld_sub)) {
      # L=1 (single effect). With L>=2 SuSiE over-fits at strong-signal loci
      # (each extra CS gets a single SNP at PIP=1 with no real LD support —
      # verified at locus_GGT_chr15_60883281: L=5 placed 5 SNPs at PIP=1
      # with pairwise r² as low as 0.002). The canonical pipeline uses
      # higher L because it has refined LD/sumstats QC; here we stay
      # conservative with L=1 for the inline fallback.
      fit <- tryCatch(
        susieR::susie_rss(z = z, R = ld_sub, L = 1,
                          coverage = 0.95, n = 343850,
                          max_iter = 500, refine = FALSE),
        error = function(e) { cat("  susie_rss failed:", e$message, "\n"); NULL })
      if (!is.null(fit)) {
        cat(sprintf("  susie_rss returned. converged=%s, max PIP=%.3f\n",
                    fit$converged, max(fit$pip)))
      }
      if (!is.null(fit) && (isTRUE(fit$converged) || max(fit$pip) >= 0.3)) {
        cs_idx <- integer(length(fit$pip))
        if (length(fit$sets$cs))
          for (k in seq_along(fit$sets$cs))
            cs_idx[fit$sets$cs[[k]]] <- k
        ukbb_pips <- data.table(
          study     = "UKBB_EUR (inline)",
          chromosome = CHR,
          position  = eur_sub$position,
          susie_pip = fit$pip,
          susie_cs  = cs_idx,
          susie_converged = TRUE)[susie_pip > 0]
        pips <- rbind(pips, ukbb_pips, fill = TRUE)
        cat(sprintf("  Inline UKBB EUR SuSiE: %d variants (max PIP=%.3f)\n",
                    nrow(ukbb_pips), max(ukbb_pips$susie_pip)))
      } else cat("  Inline SuSiE did not converge.\n")
    } else cat("  z-score length / LD dim mismatch — skipping.\n")
  } else cat("  Too few shared variants for inline SuSiE.\n")
}

# ── eQTL track data ───────────────────────────────────────────────────────────
eqtl_file <- file.path(EQTL_DIR, paste0("chr", CHR, "_marginal_summary_results.tsv"))
eqtl_track <- NULL
best_gene  <- NA_character_
if (file.exists(eqtl_file)) {
  eqtl <- fread(eqtl_file,
                select = c("Variant","CHR","POS","Beta","SE","PVAL","GeneSymbol"))
  eqtl_loc <- eqtl[CHR == CHR & POS >= WIN_START & POS <= WIN_END]
  if (nrow(eqtl_loc) > 0) {
    if (!is.na(FORCE_EQTL_GENE) && FORCE_EQTL_GENE %in% eqtl_loc$GeneSymbol) {
      best_gene <- FORCE_EQTL_GENE
    } else {
      near_lead <- eqtl_loc[abs(POS - EUR_LEAD_POS) <= 75000]
      if (!nrow(near_lead)) near_lead <- eqtl_loc
      best_gene <- near_lead[which.min(PVAL), GeneSymbol]
    }
    eqtl_track <- eqtl_loc[GeneSymbol == best_gene]
    eqtl_track[, p := as.numeric(PVAL)]
    cat(sprintf("  eQTL track: %s (%d variants)\n", best_gene, nrow(eqtl_track)))
  }
}

# ── SuSiE-COLOC top variant for the focal gene + GWAS ───────────────────────
# Resolves the variant chosen by SuSiE-COLOC as the highest joint posterior
# across (GWAS-CS × eQTL-CS) pairs — i.e. the variant most likely to be the
# shared causal between GWAS and eQTL signals. Drawn as a diamond on the
# eQTL track (and EUR GWAS track if it doesn't coincide with the lead SNP).
coloc_top_pos <- NA_integer_
coloc_top_pp  <- NA_real_
coloc_pp4     <- NA_real_
# BASE points at GWAS/finemapping (see line 42).
candidate_dirs <- c(
  file.path(BASE, "results/susie_coloc_polyfun"),
  file.path(BASE, "results/susie_coloc_1kg"),
  file.path(BASE, "results/susie_coloc"))
candidate_studies <- c(
  if (TRAIT_PAIR %in% c("ALT","AST","GGT"))
       c(paste0("UKBB_", TRAIT_PAIR), paste0("BBJ_", TRAIT_PAIR)),
  TRAIT_PAIR)  # fallback: trait-name doubles as study
for (d in candidate_dirs) {
  for (s in candidate_studies) {
    f <- file.path(d, s, sprintf("susie_coloc_chr%d.csv", CHR))
    if (!file.exists(f)) next
    tt <- fread(f)
    hit <- tt[gene == best_gene & method == "susie" &
              !is.na(PP.H4.susie) & PP.H4.susie >= 0.5]
    if (!nrow(hit)) next
    # Best CS-pair top variant by joint posterior
    setorder(hit, -top_snp_PP)
    pos <- as.integer(sub("^[0-9]+:", "", hit$top_snp[1]))
    if (!is.na(pos)) {
      coloc_top_pos <- pos
      coloc_top_pp  <- hit$top_snp_PP[1]
      coloc_pp4     <- hit$PP.H4.susie[1]
      cat(sprintf("  COLOC top SNP (%s × %s, %s): chr%d:%d  PP4=%.3f  topPP=%.3f\n",
                  best_gene, s, basename(d), CHR, coloc_top_pos,
                  coloc_pp4, coloc_top_pp))
      break
    }
  }
  if (!is.na(coloc_top_pos)) break
}

# ── Atlas: dream_logFC + dream_padj for continuous gene-track coloring ───────
# We color genes by signed bulk dream logFC modulated by significance.
#   logFC ≥ 0.05 magnitude AND padj < 0.05 → diverging blue↔red gradient
#                                            (saturated at ±LFC_SAT)
#   otherwise (n.s. or no atlas data)      → grey ("not significantly changed")
# This replaces the prior binary up/down classification; intensity of the
# colour now encodes effect size rather than a hard threshold.
LFC_SAT  <- 0.8  # saturation point for the color scale (95th-percentile shoulder)
LFC_MIN  <- 0.05 # below this magnitude, treat as essentially no change

gene_atlas <- if (file.exists(ATLAS_FILE)) {
  fread(ATLAS_FILE, select = c("human_symbol", "dream_logFC", "dream_padj"))
} else data.table(human_symbol = character(),
                  dream_logFC  = double(),
                  dream_padj   = double())

# Diverging palette: deep blue (down) → white → deep red (up)
deg_palette <- colorRampPalette(c("#1565C0", "#90CAF9", "#FFFFFF",
                                   "#EF9A9A", "#C62828"))(101)
lfc_to_color <- function(lfc, padj) {
  if (is.na(lfc) || is.na(padj) || padj >= 0.05 || abs(lfc) < LFC_MIN) {
    return("#BDBDBD")  # n.s. / no data
  }
  v <- pmin(pmax(lfc / LFC_SAT, -1), 1)
  idx <- round((v + 1) * 50) + 1   # map [-1,1] → 1..101
  deg_palette[idx]
}

# ─────────────────────────────────────────────────────────────────────────────
# Render with plotgardener
# ─────────────────────────────────────────────────────────────────────────────
out_path <- if (!is.na(OUT_OVERRIDE)) {
  dir.create(dirname(OUT_OVERRIDE), recursive = TRUE, showWarnings = FALSE)
  OUT_OVERRIDE
} else {
  file.path(OUT_DIR, paste0(LOCUS_ID, "_ld_zoom_pg.pdf"))
}
cat("Output:", out_path, "\n")

# Page geometry (inches)
PAGE_W <- 7.8
PAGE_H <- 8.4   # 2026-05-04: LD triangle removed; separate COLOC panel removed
MARGIN_L <- 0.95   # leaves room for y-axis labels
MARGIN_R <- 1.35   # leaves room for legends + per-track right-margin labels
PLOT_X <- MARGIN_L
PLOT_W <- PAGE_W - MARGIN_L - MARGIN_R

pdf(out_path, width = PAGE_W, height = PAGE_H)
pageCreate(width = PAGE_W, height = PAGE_H,
           default.units = "inches", showGuides = FALSE,
           xgrid = 0, ygrid = 0)

# Note: passing chrom/chromstart/chromend/x/width explicitly to each track instead
# of via pgParams. plotManhattan with `params = pgParams(...)` silently breaks
# vector `fill` (LD-coloring), see plotgardener internals — parseColors path
# diverges. Explicit args sidestep the issue.

# Convert "inches from top of page" → "inches from bottom" (grid native).
# plotgardener tracks use top-down inches; bare grid primitives use bottom-up.
top_to_grid <- function(y_top) PAGE_H - y_top

# LD palette for EUR Manhattan — 8-step gradient (grey → blue → green → yellow
# → orange → red), set by user 2026-05-04. Lowest r² = grey, highest = deep red.
ld_palette <- colorRampPalette(c(
  "#999999",  # 0.0  grey
  "#0E5579",  # ~0.15 dark blue
  "#0073B2",  # ~0.30 medium blue
  "#009E73",  # ~0.45 green
  "#9AD093",  # ~0.60 light green
  "#F0E442",  # ~0.75 yellow
  "#F19F3F",  # ~0.90 orange
  "#B82725"   # 1.0  deep red
))

track_label <- function(label, x = PAGE_W - MARGIN_R + 0.05, y) {
  plotText(label = label, x = x, y = y,
           fontsize = 8, fontface = "bold",
           just = c("left", "top"), default.units = "inches")
}

axis_label <- function(label, y, h) {
  plotText(label = label, rot = 90,
           x = MARGIN_L - 0.55, y = y + h/2,
           fontsize = 7, just = "center",
           default.units = "inches")
}

# Draw left + bottom axis spines for a track at (y_top, h) inches-from-top.
draw_axis_spines <- function(y, h) {
  y_top_grid    <- top_to_grid(y)
  y_bottom_grid <- top_to_grid(y + h)
  # Left (y) spine
  grid.segments(
    x0 = unit(PLOT_X, "inches"), x1 = unit(PLOT_X, "inches"),
    y0 = unit(y_bottom_grid, "inches"),
    y1 = unit(y_top_grid,    "inches"),
    gp = gpar(col = "black", lwd = 0.5))
  # Bottom (x) spine
  grid.segments(
    x0 = unit(PLOT_X,           "inches"),
    x1 = unit(PLOT_X + PLOT_W,  "inches"),
    y0 = unit(y_bottom_grid,    "inches"),
    y1 = unit(y_bottom_grid,    "inches"),
    gp = gpar(col = "black", lwd = 0.5))
}

# Track 1: EUR GWAS, LD-colored ───────────────────────────────────────────────
y1 <- 0.45; h1 <- 1.45
eur_pg <- eur_locus[!is.na(p),
  .(chrom = paste0("chr", chromosome), pos = position, p = p, r2 = r2)]
# Sort low → high r² so high-LD points draw on top.
setorder(eur_pg, r2)
ld_pal_vec   <- ld_palette(100)
eur_pg[, fill_col := ld_pal_vec[pmin(100, pmax(1, ceiling(r2 * 100)))]]
eur_pg_df  <- as.data.frame(eur_pg)
# plotgardener's plotManhattan + vector `fill` has a bug where the per-row
# colors don't reach the rendered points — they all paint with the default
# col. Workaround: emit an empty Manhattan (one color, alpha=0) just to set
# up the y-axis / coordinate system and significance line, then overlay the
# LD-coloured points ourselves with grid.points.
eur_y_max <- ceiling(max(-log10(eur_pg_df$p), na.rm = TRUE)) + 1
mh_eur <- plotManhattan(
  data = eur_pg_df,
  chrom = paste0("chr", CHR), chromstart = WIN_START, chromend = WIN_END,
  assembly = "hg19",
  fill = NA, pch = 19, cex = 0.001,    # invisible — frame only
  sigVal = 5e-8, sigLine = TRUE, sigCol = "black",
  range = c(0, eur_y_max),
  x = PLOT_X, width = PLOT_W,
  y = y1, height = h1, default.units = "inches"
)
# Manual point overlay: convert genomic pos → inch x; -log10p → inch y.
# Grid uses bottom-up coordinates while plotgardener tracks are placed
# top-down via pageCreate(); convert by subtracting from PAGE_H.
mh_eur_pos_to_x <- function(pos) {
  PLOT_X + (pos - WIN_START) / (WIN_END - WIN_START) * PLOT_W
}
mh_eur_p_to_y <- function(p_value) {
  top    <- PAGE_H - y1                    # top edge in grid (bottom-up)
  bottom <- PAGE_H - (y1 + h1)             # bottom edge in grid
  bottom + (-log10(p_value)) / eur_y_max * (top - bottom)
}
grid.points(
  x = unit(mh_eur_pos_to_x(eur_pg_df$pos), "inches"),
  y = unit(mh_eur_p_to_y(eur_pg_df$p),    "inches"),
  pch  = 19,
  size = unit(0.07, "inches"),
  gp   = gpar(col = eur_pg_df$fill_col)
)

# Highlight lead SNP with a diamond marker
lead_x <- mh_eur_pos_to_x(EUR_LEAD_POS)
lead_p <- eur_pg_df$p[which.min(abs(eur_pg_df$pos - EUR_LEAD_POS))]
lead_y <- mh_eur_p_to_y(lead_p)
grid.points(
  x = unit(lead_x, "inches"), y = unit(lead_y, "inches"),
  pch = 23, size = unit(0.13, "inches"),
  gp = gpar(col = "black", fill = "#C2185B", lwd = 0.8)
)

annoYaxis(plot = mh_eur, at = pretty(c(0, max(-log10(eur_pg$p)))), fontsize = 7)
draw_axis_spines(y1, h1)
axis_label("-log10(p)", y1, h1)
track_label("EUR GWAS", y = y1 + 0.02)
# Shared LD legend (5-bin r² + lead) — placed in the right margin centered
# vertically on the seam between EUR and EAS Manhattans so it visually serves
# both. (Drawn after Track 2 is defined; see block below the EAS track.)

# Track 2: EAS GWAS, LD-colored to EAS lead (mirrors EUR track) ───────────────
y2 <- y1 + h1 + 0.18; h2 <- 1.10
eas_pg <- eas_locus[!is.na(p),
  .(chrom = paste0("chr", chromosome), pos = position, p = p, r2 = r2)]
if (nrow(eas_pg) > 0) {
  setorder(eas_pg, r2)  # high-LD draws on top
  eas_pg[, fill_col := ld_pal_vec[pmin(100, pmax(1, ceiling(r2 * 100)))]]
  eas_pg_df <- as.data.frame(eas_pg)
  eas_y_max <- ceiling(max(-log10(eas_pg_df$p), na.rm = TRUE)) + 1
  mh_eas <- plotManhattan(
    data = eas_pg_df,
    chrom = paste0("chr", CHR), chromstart = WIN_START, chromend = WIN_END,
    assembly = "hg19",
    fill = NA, pch = 19, cex = 0.001,    # invisible — frame only
    sigVal = 5e-8, sigLine = TRUE, sigCol = "black",
    range = c(0, eas_y_max),
    x = PLOT_X, width = PLOT_W,
    y = y2, height = h2, default.units = "inches"
  )
  mh_eas_pos_to_x <- function(pos) {
    PLOT_X + (pos - WIN_START) / (WIN_END - WIN_START) * PLOT_W
  }
  mh_eas_p_to_y <- function(p_value) {
    top    <- PAGE_H - y2
    bottom <- PAGE_H - (y2 + h2)
    bottom + (-log10(p_value)) / eas_y_max * (top - bottom)
  }
  grid.points(
    x = unit(mh_eas_pos_to_x(eas_pg_df$pos), "inches"),
    y = unit(mh_eas_p_to_y(eas_pg_df$p),    "inches"),
    pch  = 19,
    size = unit(0.07, "inches"),
    gp   = gpar(col = eas_pg_df$fill_col)
  )
  # Highlight EAS lead SNP with a diamond
  eas_lead_x <- mh_eas_pos_to_x(EAS_LEAD_POS)
  eas_lead_p <- eas_pg_df$p[which.min(abs(eas_pg_df$pos - EAS_LEAD_POS))]
  eas_lead_y <- mh_eas_p_to_y(eas_lead_p)
  grid.points(
    x = unit(eas_lead_x, "inches"), y = unit(eas_lead_y, "inches"),
    pch = 23, size = unit(0.13, "inches"),
    gp = gpar(col = "black", fill = "#C2185B", lwd = 0.8)
  )
  annoYaxis(plot = mh_eas, at = pretty(c(0, max(-log10(eas_pg$p)))), fontsize = 7)
  draw_axis_spines(y2, h2)
} else {
  plotText(label = "(no EAS data in window)",
           x = PLOT_X + PLOT_W/2, y = y2 + h2/2,
           fontsize = 8, fontcolor = "grey50",
           default.units = "inches")
}
axis_label("-log10(p)", y2, h2)
track_label("EAS GWAS", y = y2 + 0.02)

# ── Shared LD legend (EUR + EAS) ─────────────────────────────────────────────
# Placed horizontally in the right margin of the EUR GWAS track. This resolves
# the awkward overlap with the "EAS GWAS" label below by keeping the legend
# compact, modern, and perfectly aligned with the Genes log2FC scale bar.
ld_x0 <- PAGE_W - MARGIN_R + 0.05
ld_y  <- y1 + 0.35
ld_w  <- 1.05

# Horizontal color scale bar (50-step gradient for smooth transition)
n_steps <- 50
for (k in 0:(n_steps - 1)) {
  step_w   <- ld_w / n_steps
  step_r2  <- k / (n_steps - 1)
  step_col <- ld_pal_vec[pmin(100, pmax(1, ceiling(step_r2 * 100)))]
  grid.rect(
    x = unit(ld_x0 + k * step_w, "inches"),
    y = unit(top_to_grid(ld_y),  "inches"),
    width  = unit(step_w + 0.005, "inches"),
    height = unit(0.10, "inches"),
    just = c("left", "center"),
    gp = gpar(col = NA, fill = step_col)
  )
}

# Scale endpoints + middle tick
plotText(label = "0.0",
         x = ld_x0, y = ld_y + 0.10,
         fontsize = 5.5, just = c("left", "top"), default.units = "inches")
plotText(label = "0.5",
         x = ld_x0 + ld_w/2, y = ld_y + 0.10,
         fontsize = 5.5, just = c("center", "top"), default.units = "inches")
plotText(label = "1.0",
         x = ld_x0 + ld_w, y = ld_y + 0.10,
         fontsize = 5.5, just = c("right", "top"), default.units = "inches")

# Scale title (above the bar)
plotText(label = "r² to lead",
         x = ld_x0 + ld_w/2, y = ld_y - 0.10,
         fontsize = 5.5, fontface = "italic", just = c("center", "bottom"),
         default.units = "inches")

# Lead SNP indicator (below the scale labels)
lead_y_offset <- ld_y + 0.30
grid.points(
  x = unit(ld_x0 + 0.05, "inches"),
  y = unit(top_to_grid(lead_y_offset), "inches"),
  pch = 23, size = unit(0.10, "inches"),
  gp = gpar(col = "black", fill = "#C2185B", lwd = 0.6)
)
plotText(label = "Lead SNP",
         x = ld_x0 + 0.15, y = lead_y_offset,
         fontsize = 5.5, just = c("left", "center"),
         default.units = "inches")

# Track 2.5 (LD triangular heatmap) dropped 2026-05-04 — user request.
# The block was previously here (~125 lines, used `y_ld` and `h_ld`). Track 3
# (PIP) now starts directly after Track 2 (EAS).
if (FALSE) {
y_ld <- y2 + h2 + 0.18; h_ld <- 1.05
cat("[ld-triangle] building rotated LD heatmap ...\n")

# Thin LD panel uniformly to ~200 variants (~20K pairs)
N_target  <- 200L
thin_step <- max(1L, floor(length(ld_pos) / N_target))
keep_idx  <- seq(1L, length(ld_pos), by = thin_step)
ld_thin   <- ld_mat[keep_idx, keep_idx, drop = FALSE]
pos_thin  <- ld_pos[keep_idx]
n_snps    <- length(pos_thin)
cat(sprintf("  thinned LD panel: %d variants → %d pairs\n",
            n_snps, n_snps * (n_snps + 1) / 2))

if (n_snps >= 2) {
  # Cream → light green → cyan → navy palette (YlGnBu)
  ld_tri_pal_fn <- colorRampPalette(c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB",
                                       "#41B6C4", "#1D91C0", "#225EA8", "#081D58"))

  win_range       <- WIN_END - WIN_START
  max_depth_shown <- 0.20                # show only top 20% near-diag wedge

  # Build a high-resolution image (NX × NY pixels) where each pixel covers a
  # (x_centre, depth) cell. For each upper-triangle pair (i, j), drop its r²
  # value into the matching pixel — the image is then drawn with grid.raster
  # which fills the wedge cleanly.
  NX <- 1600L
  NY <- 420L
  img <- matrix(NA_real_, nrow = NY, ncol = NX)

  pairs_dt <- CJ(i = seq_len(n_snps), j = seq_len(n_snps))[i <= j]
  pairs_dt[, `:=`(p_i = pos_thin[i], p_j = pos_thin[j])]
  pairs_dt[, r2_pair    := pmin(1, ld_thin[cbind(i, j)]^2)]
  pairs_dt[, x_centre   := (p_i + p_j) / 2]
  pairs_dt[, depth_frac := (p_j - p_i) / win_range]
  pairs_dt <- pairs_dt[depth_frac <= max_depth_shown & !is.na(r2_pair)]

  pairs_dt[, ix := pmin(NX, pmax(1L,
              ceiling((x_centre - WIN_START) / win_range * NX)))]
  pairs_dt[, iy := pmin(NY, pmax(1L,
              ceiling(depth_frac / max_depth_shown * NY)))]

  # Splat each pair into a footprint of pixels so neighbouring pairs tile the
  # wedge without gaps. The original ggplot version used ~8× the average
  # spacing (via `max(thin_step, 8)`); we use 5× here to get a clean
  # continuous fill while still preserving block boundaries (pmax keeps the
  # high-LD blocks visible against neighbouring low-r² pairs).
  fp_x_px <- max(4L, ceiling(NX * (mean(diff(pos_thin)) / win_range) * 5.0))
  fp_y_px <- max(4L, ceiling(NY * (mean(diff(pos_thin)) / (max_depth_shown * win_range)) * 5.0))

  for (k in seq_len(nrow(pairs_dt))) {
    cx <- pairs_dt$ix[k]; cy <- pairs_dt$iy[k]; rv <- pairs_dt$r2_pair[k]
    xs <- max(1L, cx - fp_x_px) : min(NX, cx + fp_x_px)
    ys <- max(1L, cy - fp_y_px) : min(NY, cy + fp_y_px)
    cur <- img[ys, xs]
    new <- pmax(cur, rv, na.rm = TRUE)
    img[ys, xs] <- new
  }

  # Map r² values → palette colours; NA pixels (outside wedge) → very light grey.
  to_color <- function(v) {
    out <- rep(NA_character_, length(v))
    finite <- !is.na(v)
    out[!finite] <- "#FAFAFA"
    if (any(finite)) {
      idx <- pmin(101L, pmax(1L, round(v[finite] * 100) + 1L))
      out[finite] <- ld_tri_pal_fn(101)[idx]
    }
    out
  }
  raster_mat <- matrix(to_color(as.numeric(img)),
                       nrow = NY, ncol = NX)

  grid.raster(
    image = raster_mat,
    x = unit(PLOT_X,                       "inches"),
    y = unit(top_to_grid(y_ld + h_ld),     "inches"),
    width  = unit(PLOT_W,                  "inches"),
    height = unit(h_ld,                    "inches"),
    just = c("left", "bottom"),
    interpolate = FALSE
  )

  # Re-use the same lookup for the legend.
  ld_tri_lut <- ld_tri_pal_fn(101)

  # Top edge tick at the lead position
  lead_x_in <- PLOT_X + (EUR_LEAD_POS - WIN_START) / win_range * PLOT_W
  grid.segments(
    x0 = unit(lead_x_in, "inches"),
    y0 = unit(top_to_grid(y_ld) + 0.02, "inches"),
    x1 = unit(lead_x_in, "inches"),
    y1 = unit(top_to_grid(y_ld) - 0.05, "inches"),
    gp = gpar(col = "#C2185B", lwd = 0.8)
  )

  # Compact LD-triangle colour scale on the right margin
  # Place the LD r² scale bar directly below the "LD (EUR r²)" track label.
  ld_tri_x0 <- PAGE_W - MARGIN_R + 0.05
  ld_tri_y  <- y_ld + 0.32
  ld_tri_w  <- 1.05
  for (k in 0:100) {
    step_w <- ld_tri_w / 101
    grid.rect(
      x = unit(ld_tri_x0 + k * step_w, "inches"),
      y = unit(top_to_grid(ld_tri_y),  "inches"),
      width  = unit(step_w + 0.005,    "inches"),
      height = unit(0.10,              "inches"),
      just = c("left", "centre"),
      gp = gpar(col = NA, fill = ld_tri_lut[k + 1])
    )
  }
  plotText(label = "0",   x = ld_tri_x0,            y = ld_tri_y + 0.10,
           fontsize = 5.5, just = c("left", "top"), default.units = "inches")
  plotText(label = "1",   x = ld_tri_x0 + ld_tri_w, y = ld_tri_y + 0.10,
           fontsize = 5.5, just = c("right", "top"), default.units = "inches")
  plotText(label = "r²",  x = ld_tri_x0 + ld_tri_w/2,
           y = ld_tri_y - 0.10,
           fontsize = 5.5, fontface = "italic",
           just = c("centre", "bottom"), default.units = "inches")
} else {
  cat("  too few LD variants for triangle — skipping.\n")
}
track_label("LD (EUR r²)", y = y_ld + 0.02)
}  # end of dropped Track 2.5 block

# Track 3: PIP lollipop (EUR + EAS overlaid) ──────────────────────────────────
y3 <- y2 + h2 + 0.20; h3 <- 1.30
# Use grid primitives to draw lollipops on top of the genomic axis. grid is
# bottom-up while plotgardener is top-down; helpers below convert.
# Coordinates of the track in inches FROM TOP:
pip_top_top    <- y3 + 0.10
pip_baseline_t <- y3 + h3 - 0.10
pip_baseline <- top_to_grid(pip_baseline_t)
pip_top      <- top_to_grid(pip_top_top)
pos_to_inches <- function(pos) {
  PLOT_X + (pos - WIN_START) / (WIN_END - WIN_START) * PLOT_W
}
pip_to_inches <- function(p) {
  # p ∈ [0,1] → y in grid (bottom-up) inches between baseline and top
  pip_baseline + p * (pip_top - pip_baseline)
}

# Frame: baseline + threshold line (PIP = 0.9) + left spine
grid.segments(
  x0 = unit(PLOT_X,            "inches"),
  x1 = unit(PLOT_X + PLOT_W,   "inches"),
  y0 = unit(pip_baseline,      "inches"),
  y1 = unit(pip_baseline,      "inches"),
  gp = gpar(col = "black", lwd = 0.6)
)
grid.segments(
  x0 = unit(PLOT_X,            "inches"),
  x1 = unit(PLOT_X,            "inches"),
  y0 = unit(pip_baseline,      "inches"),
  y1 = unit(pip_top,           "inches"),
  gp = gpar(col = "black", lwd = 0.6)
)
grid.segments(
  x0 = unit(PLOT_X,                          "inches"),
  x1 = unit(PLOT_X + PLOT_W,                 "inches"),
  y0 = unit(pip_to_inches(0.9),              "inches"),
  y1 = unit(pip_to_inches(0.9),              "inches"),
  gp = gpar(col = "grey50", lty = 2, lwd = 0.4)
)

# y-axis ticks at 0, 0.5, 1.0
for (val in c(0, 0.5, 1.0)) {
  grid.segments(
    x0 = unit(PLOT_X - 0.05,        "inches"),
    x1 = unit(PLOT_X,               "inches"),
    y0 = unit(pip_to_inches(val),   "inches"),
    y1 = unit(pip_to_inches(val),   "inches"),
    gp = gpar(col = "black", lwd = 0.5)
  )
  grid.text(
    label = format(val, nsmall = 1),
    x = unit(PLOT_X - 0.08,         "inches"),
    y = unit(pip_to_inches(val),    "inches"),
    just = "right",
    gp = gpar(fontsize = 6)
  )
}

# Lollipops per ancestry. UKBB_* = blue (EUR); BBJ_* = orange (EAS); other
# studies = blue (treat as EUR by default).
study_color <- function(s) {
  if (grepl("^BBJ_", s))     return("#D55E00")
  if (grepl("^UKBB", s))     return("#1565C0")
  if (grepl("_EAS$", s))     return("#D55E00")
  if (grepl("EUR.*inline", s)) return("#1565C0")
  "#1565C0"
}
draw_susie_pip <- function(dt) {
  if (!nrow(dt)) return(invisible())
  for (i in seq_len(nrow(dt))) {
    xi   <- pos_to_inches(dt$position[i])
    pv   <- dt$susie_pip[i]
    yi   <- pip_to_inches(pv)
    col  <- study_color(dt$study[i])
    in_cs <- !is.na(dt$susie_cs[i]) && dt$susie_cs[i] > 0
    grid.segments(x0 = unit(xi, "inches"), x1 = unit(xi, "inches"),
                  y0 = unit(pip_baseline, "inches"),
                  y1 = unit(yi, "inches"),
                  gp = gpar(col = "grey85", lwd = 0.3))
    grid.points(x = unit(xi, "inches"), y = unit(yi, "inches"),
                pch = if (in_cs) 19 else 21,
                size = unit(if (in_cs) 0.10 else 0.05, "inches"),
                gp = gpar(col = col, fill = if (in_cs) col else "white",
                          lwd = 0.5))
  }
}
draw_susie_pip(pips)

axis_label("PIP", y3, h3)
track_label("SuSiE finemap PIP", y = y3 + 0.02)

# Inline mini-legend
pip_legend_x <- PLOT_X + 0.10
pip_legend_y <- y3 + 0.10
present_studies <- unique(pips$study)
ukbb_present <- any(grepl("^UKBB", present_studies))
bbj_present  <- any(grepl("^BBJ_", present_studies))
xc <- pip_legend_x
if (ukbb_present) {
  grid.points(x = unit(xc, "inches"),
              y = unit(top_to_grid(pip_legend_y), "inches"),
              pch = 19, size = unit(0.07, "inches"), gp = gpar(col = "#1565C0"))
  plotText(label = "UKBB EUR (in CS)", x = xc + 0.06, y = pip_legend_y,
           fontsize = 5.5, fontcolor = "grey25",
           just = c("left", "center"), default.units = "inches")
  xc <- xc + 0.95
}
if (bbj_present) {
  grid.points(x = unit(xc, "inches"),
              y = unit(top_to_grid(pip_legend_y), "inches"),
              pch = 19, size = unit(0.07, "inches"), gp = gpar(col = "#D55E00"))
  plotText(label = "BBJ EAS (in CS)", x = xc + 0.06, y = pip_legend_y,
           fontsize = 5.5, fontcolor = "grey25",
           just = c("left", "center"), default.units = "inches")
  xc <- xc + 0.92
}
grid.points(x = unit(xc, "inches"),
            y = unit(top_to_grid(pip_legend_y), "inches"),
            pch = 21, size = unit(0.05, "inches"),
            gp = gpar(col = "grey50", fill = "white", lwd = 0.5))
plotText(label = "not in CS", x = xc + 0.06, y = pip_legend_y,
         fontsize = 5.5, fontcolor = "grey25",
         just = c("left", "center"), default.units = "inches")

# (Separate SuSiE-COLOC PP.H4 panel removed 2026-05-04 — diamond restored on
#  the eQTL track below at the COLOC top variant.)

# Track 4: eQTL Manhattan (focal gene) ────────────────────────────────────────
y4 <- y3 + h3 + 0.18; h4 <- 0.95
if (!is.null(eqtl_track) && nrow(eqtl_track) > 0) {
  eqtl_pg <- eqtl_track[!is.na(p) & p > 0,
    .(chrom = paste0("chr", CHR), pos = POS, p = p)]
  mh_eqtl <- plotManhattan(
    data = eqtl_pg,
    chrom = paste0("chr", CHR), chromstart = WIN_START, chromend = WIN_END,
    assembly = "hg19",
    fill = "#E07B39",
    pch = 19, cex = 0.40,
    sigVal = 1e-5, sigLine = FALSE,
    range = c(0, ceiling(max(-log10(eqtl_pg$p), na.rm = TRUE)) + 1),
    x = PLOT_X, width = PLOT_W,
    y = y4, height = h4, default.units = "inches"
  )
  annoYaxis(plot = mh_eqtl, at = pretty(c(0, max(-log10(eqtl_pg$p)))), fontsize = 7)
  draw_axis_spines(y4, h4)
  axis_label("-log10(p)", y4, h4)
  track_label(sprintf("eQTL: %s", best_gene), y = y4 + 0.02)

  # Diamond at SuSiE-COLOC top variant (highest joint posterior).
  if (!is.na(coloc_top_pos) &&
      coloc_top_pos >= WIN_START && coloc_top_pos <= WIN_END) {
    eqtl_y_max <- ceiling(max(-log10(eqtl_pg$p), na.rm = TRUE)) + 1
    coloc_p <- eqtl_pg$p[which.min(abs(eqtl_pg$pos - coloc_top_pos))]
    coloc_x_in <- PLOT_X + (coloc_top_pos - WIN_START) /
                  (WIN_END - WIN_START) * PLOT_W
    coloc_y_top <- PAGE_H - y4
    coloc_y_bot <- PAGE_H - (y4 + h4)
    coloc_y_in  <- coloc_y_bot + (-log10(coloc_p)) /
                   eqtl_y_max * (coloc_y_top - coloc_y_bot)
    grid.points(
      x = unit(coloc_x_in, "inches"), y = unit(coloc_y_in, "inches"),
      pch = 23, size = unit(0.13, "inches"),
      gp = gpar(col = "black", fill = "#FFD600", lwd = 0.8))
    plotText(label = sprintf("COLOC top SNP\nPP4 = %.2f", coloc_pp4),
             x = coloc_x_in + 0.10, y = y4 + 0.10,
             fontsize = 6, fontface = "italic",
             just = c("left", "top"), default.units = "inches",
             fontcolor = "grey20")
  }
}

# Track 5: gene track (plotgardener TxDb), colored by bulk dream logFC ────────
y5 <- y4 + h4 + 0.20; h5 <- 1.00

# Pre-fetch the genes plotGenes will draw, so we can hand each a logFC color.
locus_gr     <- GRanges(seqnames = paste0("chr", CHR),
                        ranges   = IRanges(WIN_START, WIN_END))
genes_in_win <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                       filter = list(tx_chrom = paste0("chr", CHR))))
genes_in_win <- subsetByOverlaps(genes_in_win, locus_gr)
gene_symbols <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = genes_in_win$gene_id,
         column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
)
gene_symbols <- unique(gene_symbols[!is.na(gene_symbols)])

# Map each gene → continuous color via dream_logFC + dream_padj.
gene_colors <- setNames(character(length(gene_symbols)), gene_symbols)
for (g in gene_symbols) {
  hit <- gene_atlas[human_symbol == g][1]
  gene_colors[g] <- if (nrow(hit) && !is.na(hit$dream_logFC)) {
    lfc_to_color(hit$dream_logFC, hit$dream_padj)
  } else "#BDBDBD"
}
gene_hl <- data.frame(gene  = names(gene_colors),
                      color = unname(gene_colors),
                      stringsAsFactors = FALSE)

plotGenes(
  chrom = paste0("chr", CHR), chromstart = WIN_START, chromend = WIN_END,
  assembly = "hg19",
  fill = c("#9E9E9E", "#9E9E9E"),
  fontcolor = c("#424242", "#424242"),
  geneHighlights = gene_hl,
  geneBackground = "grey85",
  fontsize = 7, strandLabels = TRUE,
  x = PLOT_X, width = PLOT_W,
  y = y5, height = h5, default.units = "inches"
)
track_label("Genes", y = y5 + 0.02)

# Inline horizontal logFC scale bar in the right margin
# Place the logFC scale bar directly below the "Genes" track label.
sb_x0 <- PAGE_W - MARGIN_R + 0.05
sb_y  <- y5 + 0.32
sb_w  <- 1.05
n_steps <- 50
for (k in 0:(n_steps - 1)) {
  step_w <- sb_w / n_steps
  step_lfc <- (-LFC_SAT) + (2 * LFC_SAT) * (k / (n_steps - 1))
  step_col <- deg_palette[round((step_lfc / LFC_SAT + 1) * 50) + 1]
  grid.rect(
    x = unit(sb_x0 + k * step_w, "inches"),
    y = unit(top_to_grid(sb_y),  "inches"),
    width  = unit(step_w + 0.005, "inches"),
    height = unit(0.10, "inches"),
    just = c("left", "center"),
    gp = gpar(col = NA, fill = step_col)
  )
}
# Scale endpoints + zero tick
plotText(label = sprintf("-%.1f", LFC_SAT),
         x = sb_x0, y = sb_y + 0.10,
         fontsize = 5.5, just = c("left","top"), default.units = "inches")
plotText(label = "0",
         x = sb_x0 + sb_w/2, y = sb_y + 0.10,
         fontsize = 5.5, just = c("center","top"), default.units = "inches")
plotText(label = sprintf("+%.1f", LFC_SAT),
         x = sb_x0 + sb_w, y = sb_y + 0.10,
         fontsize = 5.5, just = c("right","top"), default.units = "inches")
plotText(label = "log2FC", x = sb_x0 + sb_w/2, y = sb_y - 0.10,
         fontsize = 5.5, fontface = "italic", just = c("center","bottom"),
         default.units = "inches")
# Caveat: n.s. genes are grey (a single neutral indicator)
grid.rect(
  x = unit(sb_x0,                 "inches"),
  y = unit(top_to_grid(sb_y + 0.30), "inches"),
  width  = unit(0.10, "inches"),
  height = unit(0.10, "inches"),
  just = c("left","center"),
  gp = gpar(col = NA, fill = "#BDBDBD")
)
plotText(label = "n.s.",
         x = sb_x0 + 0.13, y = sb_y + 0.30,
         fontsize = 5.5, just = c("left","center"),
         default.units = "inches")

# Track 6: genome label (Mb scale) ────────────────────────────────────────────
y6 <- y5 + h5 + 0.05
plotGenomeLabel(
  chrom = paste0("chr", CHR), chromstart = WIN_START, chromend = WIN_END,
  assembly = "hg19",
  scale = "Mb", commas = TRUE, sequence = FALSE,
  fontsize = 8,
  x = PLOT_X, y = y6, length = PLOT_W, default.units = "inches"
)

pageGuideHide()
dev.off()

cat("\nDone:", out_path, "\n")
