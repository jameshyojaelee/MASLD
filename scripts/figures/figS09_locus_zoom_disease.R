#!/usr/bin/env Rscript
# figS09_locus_zoom_disease.R — single-ancestry (EUR) locus zoom for colocalizations
# driven by a DISEASE-ENDPOINT (NAFLD/NASH) or IMAGING (PDFF/liver-fat) GWAS,
# rather than the liver-enzyme proxies (ALT/AST/GGT) handled by
# figS09_locus_zoom_pg.R.
#
# Why a separate script: figS09_locus_zoom_pg.R is hardwired to the EUR(UKBB-enzyme)
# + EAS(BBJ-enzyme) SuSiEx shared-locus paradigm (EUR_SS = UKBB_<trait>, locus must
# exist in shared_loci.csv). The disease/PDFF colocs are EUR-only and absent from
# shared_loci.csv, so they need a single-ancestry renderer. This reuses figS09's
# verified plotgardener coordinate math, LD loader, COLOC-diamond logic, and the
# (user-approved 2026-06-19) custom transcript gene track verbatim.
#
# Tracks (top -> bottom):
#   1. EUR disease/PDFF GWAS Manhattan, LD-colored (r^2 to lead)
#   2. SuSiE finemap PIP lollipops (EUR; stored in combined_finemapping_1kg.csv)
#   3. eQTL Manhattan (focal gene) + gold SuSiE-COLOC top-variant diamond
#   4. gene track (TxDb hg19, exon/intron transcript model, DEG-logFC colored)
#   5. plotGenomeLabel ruler (Mb)
#
# Usage:  Rscript figS09_locus_zoom_disease.R <gene> [out_override]
#   <gene> in {IL18R1, SHMT1, C2orf16, F13B, ERCC2} (config table below)
# Run inside the `plotgardener` micromamba env (R 4.4.3 + plotgardener 1.12.0).
#
# Output: figures/main/fig2_genetics/locus_zoom/<TRAIT>_<GENE>.pdf (2026-07-01;
# was fig2_genetics/panels/locuszoom_disease_<GENE>.pdf). Layout compacted
# (~5.4 x 4.2in, was 7.8 x 8.7in) with a uniform 6pt font floor.

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

# ── plotgardener: plain chromosome label on the genome ruler ──────────────────
# plotgardener hardcodes `object$gp$fontface <- "bold"` on the chrN label inside
# the non-exported plotChromGenomeLabel(); no argument of plotGenomeLabel() can
# reach it (setGP's fontface is overwritten downstream). Project style is 6 pt
# plain for every text element, so patch that one assignment. Face only — no
# position, size, or content changes.
local({
  ns <- asNamespace("plotgardener")
  src <- deparse(get("plotChromGenomeLabel", envir = ns))
  n_bold <- sum(grepl('object$gp$fontface <- "bold"', src, fixed = TRUE))
  stopifnot(n_bold > 0L)
  src <- gsub('object$gp$fontface <- "bold"', 'object$gp$fontface <- "plain"',
              src, fixed = TRUE)
  f <- eval(parse(text = paste(src, collapse = "\n")))
  environment(f) <- ns
  assignInNamespace("plotChromGenomeLabel", f, ns = "plotgardener")
})

# ── Per-gene config ───────────────────────────────────────────────────────────
# Each entry: the colocalizing EUR disease/PDFF GWAS, its hg19 sumstats file, the
# coloc lead (hg19), the focal eQTL gene, and a human-readable trait label. All
# verified against susie_coloc_all_gwas.csv + on-disk sumstats (build = hg19).
#
# Figure roster (2026-06-22): KEEP ERCC2, F13B (PDFF/imaging) + IL18R1 (the disease-
# DIAGNOSIS example: FinnGen NAFLD, noncoding regulatory IL18R1 effector — RIGHT
# mechanism class; caveat = suggestive/sub-genome-wide GWAS + coloc cannot resolve
# the single shared variant, SNP.PP.H4 0.033, so its gold lollipop sits near 0).
# DROPPED C2orf16 (config kept for provenance/repro): lead 2:27,730,940 is the GCKR
# P446L *coding* missense (rs1260326); coloc assigns the signal to C2orf16 only
# because GCKR's protein change is invisible to eQTL coloc (GCKR mRNA-eQTL PP4 0.009)
# — a coding-class locus (GCKR, like PNPLA3/TM6SF2) miscast as an expression effector.
CONFIG <- list(
  IL18R1 = list(
    eur_ss = "FinnGen_NAFLD_reformatted_hg19.tsv",       study = "FinnGen_NAFLD",
    chr = 2L,  lead = 102979658L, eqtl = "IL18R1", n_gwas = 400000L,
    gwas_label = "FinnGen NAFLD diagnosis (EUR)", short_lab = "NAFLD GWAS", trait = "NAFLD diagnosis"),
  # SHMT1: the only NASH-diagnosis coloc (FinnGen NASH). Complements IL18R1 — its
  # COLOC resolves the shared variant SHARPLY (SNP.PP.H4 0.951, gold lollipop near
  # the top) where IL18R1's is near 0, and the SHMT1 liver eQTL is razor-sharp
  # (-log10p 64.6). CAVEATS (carried in the source CSV / caption): (1) FinnGen NASH
  # signal is sub-suggestive in-window (peak -log10p 4.83, below 1e-5) -> GWAS track
  # renders flat; (2) ABF dissents from SuSiE (PP.H4.abf 0.012 / PP.H3.abf 0.755 vs
  # PP.H4.susie 0.879) — expected at this multi-signal 17p11.2 segmental-duplication
  # locus where ABF's single-causal-variant model is misspecified; (3) SHMT1 is a
  # human non-DEG (mouse DEG, LFC -0.74 padj 0.008) -> gene track renders gray.
  SHMT1 = list(
    eur_ss = "FinnGen_NASH_reformatted_hg19.tsv",        study = "FinnGen_NASH",
    chr = 17L, lead = 18273850L, eqtl = "SHMT1", n_gwas = 400000L,
    gwas_label = "FinnGen NASH diagnosis (EUR)", short_lab = "NASH GWAS", trait = "NASH diagnosis"),
  C2orf16 = list(
    eur_ss = "2020_32298765_NAFLD_EUR_preprocessed.tsv", study = "2020_32298765_NAFLD_EUR",
    chr = 2L,  lead = 27730940L,  eqtl = "C2orf16", n_gwas = 778614L,
    gwas_label = "NAFLD diagnosis GWAS (EUR, 2020)", short_lab = "NAFLD GWAS", trait = "NAFLD diagnosis"),
  F13B = list(
    eur_ss = "2021_34128465_PDFF_EUR_preprocessed.tsv",  study = "2021_34128465_PDFF_EUR",
    chr = 1L,  lead = 196973183L, eqtl = "F13B", n_gwas = 32976L,
    gwas_label = "Liver fat / PDFF GWAS (EUR, 2021)", short_lab = "PDFF GWAS", trait = "MRI liver fat (PDFF)"),
  ERCC2 = list(
    eur_ss = "2021_34128465_PDFF_EUR_preprocessed.tsv",  study = "2021_34128465_PDFF_EUR",
    chr = 19L, lead = 45855262L,  eqtl = "ERCC2", n_gwas = 32976L,
    gwas_label = "Liver fat / PDFF GWAS (EUR, 2021)", short_lab = "PDFF GWAS", trait = "MRI liver fat (PDFF)"),
  # RORA (2026-07-06): MVP R4 NAFLD-diagnosis EUR coloc — a NONCODING regulatory
  # effector at 15q22 (RORA nuclear-receptor). COLOC (susie_coloc/MVP_NAFLD_EUR,
  # chr15) PP.H4.susie=0.983 / PP.H4.abf=0.983, shared variant 15:60,878,030
  # (SNP.PP.H4 0.375); replicated in 2023_36280732_NAFLD_UKBB_EUR (abf 0.772).
  # GWAS lead 15:60,899,031 (p=7.6e-8, just below the 5e-8 line). Explicit window
  # (win_start/win_end) matches the 15q22 credible-set span 60.63-61.13 Mb.
  RORA = list(
    eur_ss = "MVP_NAFLD_EUR_reformatted_hg19.tsv",       study = "MVP_NAFLD_EUR",
    chr = 15L, lead = 60899031L, eqtl = "RORA", n_gwas = 216000L,
    win_start = 60630000L, win_end = 61130000L,
    gwas_label = "MVP NAFLD diagnosis (EUR)", short_lab = "NAFLD GWAS", trait = "NAFLD diagnosis"),
  # FABP1 (2026-07-07): UKBB ALT (liver enzyme) EUR coloc — EUR-only effector at 2p11.
  # gene_level_coloc PP.H4.susie 0.998 (UKBB_ALT); Broadaway FABP1 cis-eQTL present on
  # chr2 (~88.24 Mb). Rendered with the same compact layout as the main locus panels,
  # cairo_pdf/Helvetica). Lead 2:88,424,066 (hg19); default +/-200 kb window.
  FABP1 = list(
    eur_ss = "UKBB_ALT_reformatted_hg19.tsv",            study = "UKBB_ALT",
    chr = 2L,  lead = 88424066L, eqtl = "FABP1", n_gwas = 400000L,
    gwas_label = "UKBB ALT (EUR)", short_lab = "ALT GWAS", trait = "ALT (liver enzyme)")
)

HALF_WIN <- 200000L   # +/- window around the lead (hg19 bp)

args <- commandArgs(trailingOnly = TRUE)
# Paired-panel label sharing (Fig2G left / Fig2H right): draw the y-axis titles on
# the LEFT panel only, the shared right-side LEGEND on the RIGHT panel only. Numeric
# ticks + per-locus titles stay on both. Set via env at render time.
NO_YTITLE <- nzchar(Sys.getenv("LZ_NO_YTITLE"))       # suppress the left-edge y-axis TITLES (-log10(p)/PIP)
NO_LEGEND <- nzchar(Sys.getenv("LZ_NO_LEGEND"))       # suppress the shared right-margin legend (r²/Lead SNP/PIP keys/log2FC)
NO_TRACKLABEL <- nzchar(Sys.getenv("LZ_NO_TRACKLABEL")) # suppress the right-edge track-type labels (GWAS/SuSiE PIP/eQTL/Genes)
# Opt-in Figure 2 G/H layout. This keeps the general-purpose locus renderer
# backward compatible while making the two main panels geometrically identical.
# The shared keys live in Figure 2I, so matched-main panels never reserve an
# internal legend or right-hand track-label strip.
MATCHED_MAIN <- nzchar(Sys.getenv("LZ_MATCHED_MAIN"))
if (MATCHED_MAIN) {
  NO_YTITLE <- FALSE
  NO_LEGEND <- TRUE
  NO_TRACKLABEL <- TRUE
}
GENE <- if (length(args) >= 1) args[1] else "IL18R1"
if (!GENE %in% names(CONFIG)) stop("Unknown gene: ", GENE, " (have: ",
                                   paste(names(CONFIG), collapse = ", "), ")")
cfg <- CONFIG[[GENE]]
OUT_OVERRIDE <- if (length(args) >= 2) args[2] else NA_character_

CHR          <- cfg$chr
LEAD_POS     <- cfg$lead
EUR_LEAD_POS <- LEAD_POS
# Window: explicit per-config win_start/win_end when supplied (e.g. RORA, to match
# the 15q22 credible-set span), else the default LEAD +/- HALF_WIN.
WIN_START    <- if (!is.null(cfg$win_start)) cfg$win_start else LEAD_POS - HALF_WIN
WIN_END      <- if (!is.null(cfg$win_end))   cfg$win_end   else LEAD_POS + HALF_WIN
EQTL_GENE    <- cfg$eqtl
STUDY        <- cfg$study

cat(sprintf("Gene: %s | GWAS: %s | chr%d:%d-%d | lead %d | eQTL %s\n",
            GENE, STUDY, CHR, WIN_START, WIN_END, LEAD_POS, EQTL_GENE))

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
PROJ <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
COLOC_INPUT <- Sys.getenv(
  "FIG2_COLOC_INPUT",
  file.path(BASE, "results/susie_coloc/susie_coloc_all_gwas.csv")
)
normalize_page <- function(path, w, h) {
  python <- Sys.getenv(
    "MASLD_FIGURE_PYTHON",
    "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
  )
  status <- system2(
    python,
    c(shQuote(file.path(PROJ, "scripts/figures/normalize_pdf_page_box.py")),
      shQuote(path), format(w, trim = TRUE), format(h, trim = TRUE))
  )
  if (!identical(status, 0L)) stop("Could not normalize PDF page box: ", path)
}
EUR_SS     <- file.path(BASE, "data/sumstats", cfg$eur_ss)
EQTL_DIR   <- file.path(PROJ, "data/broadaway_eqtl")
ATLAS_FILE <- file.path(PROJ, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
PANEL_DIR  <- file.path(PROJ, "figures/main/fig2_genetics/locus_zoom")
stopifnot(file.exists(EUR_SS))

source(file.path(BASE, "src/finemapping_functions.R"))
LD_BASE   <- get_ld_base_dir("EUR")
LD_BLOCKS <- file.path(LD_BASE, "approx_LD_blocks.txt")

TRAIT_TOKEN <- sub(" GWAS", "", cfg$short_lab)   # e.g. "NAFLD" — matches figS09_locus_zoom_pg.R's <TRAIT>_<GENE> scheme
out_path <- if (!is.na(OUT_OVERRIDE)) OUT_OVERRIDE else
  file.path(PANEL_DIR, paste0(TRAIT_TOKEN, "_", GENE, ".pdf"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

# ── EUR sumstats ──────────────────────────────────────────────────────────────
cat("Loading EUR sumstats...\n")
eur_ss <- fread(EUR_SS, select = c("chromosome","position","allele1","allele2","beta","se","pval"))
eur_locus <- eur_ss[chromosome == CHR & position >= WIN_START & position <= WIN_END]
eur_locus[, p := as.numeric(pval)]
rm(eur_ss); gc(verbose = FALSE)
cat("  EUR locus variants:", nrow(eur_locus), "\n")

# ── SuSiE PIPs (stored, single-ancestry EUR; multi-source) ───────────────────
# No single finemapping aggregate covers every disease/PDFF study (the canonical
# coloc came from the polyfun arm). Try the 1kg aggregate, then the coloc-targeted
# aggregate; if both are empty for this study/window, fall back to inline SuSiE
# (below, after the LD matrix is loaded) — exactly as figS09_locus_zoom_pg.R does
# for UKBB EUR. `pip_source` records which path populated the track.
load_stored_pips <- function(fname) {
  f <- file.path(BASE, "results", fname)
  if (!file.exists(f)) return(data.table())
  cols <- names(fread(f, nrows = 1))
  fm <- fread(f, select = intersect(c("study","chromosome","position","susie_pip",
                                       "susie_cs","susie_converged"), cols))
  fm[study == STUDY & chromosome == CHR & position >= WIN_START & position <= WIN_END &
     susie_converged == TRUE & !is.na(susie_pip) & susie_pip > 0]
}
pips <- load_stored_pips("combined_finemapping_1kg.csv")
pip_source <- "1kg"
if (!nrow(pips)) { pips <- load_stored_pips("combined_finemapping_coloc_targeted.csv"); pip_source <- "coloc_targeted" }
cat(sprintf("  SuSiE finemap PIPs (stored %s): %d variants (study %s)\n", pip_source, nrow(pips), STUDY))

# ── LD r^2 to lead (EUR 1kg panel) ───────────────────────────────────────────
load_ld_for_window <- function(ld_base, ld_blocks_file) {
  blocks_tbl <- fread(ld_blocks_file)
  matching   <- blocks_tbl[chr == CHR & start < WIN_END & stop > WIN_START]
  ld_chunks <- lapply(seq_len(nrow(matching)), function(k) {
    blk <- file.path(ld_base, paste0("chr", CHR),
                     paste0(matching$start[k], ".", matching$stop[k]),
                     paste0(matching$start[k], ".", matching$stop[k]))
    bim_file <- paste0(blk, ".bim"); ld_file <- paste0(blk, ".ld")
    if (!file.exists(bim_file) || !file.exists(ld_file)) return(NULL)
    bim <- fread(bim_file, col.names = c("chr","rsid","dk","pos","alt","ref"))
    idx <- which(bim$pos >= WIN_START & bim$pos <= WIN_END)
    if (!length(idx)) return(NULL)
    ld_full <- as.matrix(fread(ld_file)); ld_sub <- ld_full[idx, idx, drop = FALSE]
    rownames(ld_sub) <- colnames(ld_sub) <- bim$pos[idx]
    rm(ld_full); gc(verbose = FALSE)
    list(ld = ld_sub, pos = bim$pos[idx])
  })
  ld_chunks <- ld_chunks[!sapply(ld_chunks, is.null)]
  if (!length(ld_chunks)) return(NULL)
  if (length(ld_chunks) == 1) list(mat = ld_chunks[[1]]$ld, pos = ld_chunks[[1]]$pos)
  else {
    positions <- unlist(lapply(ld_chunks, `[[`, "pos"))
    mat <- as.matrix(Matrix::bdiag(lapply(ld_chunks, `[[`, "ld")))
    rownames(mat) <- colnames(mat) <- positions
    list(mat = mat, pos = positions)
  }
}
r2_to_lead <- function(ld_obj, lead_pos) {
  if (is.null(ld_obj)) return(data.table(position = integer(0), r2 = double(0)))
  lead_row <- which.min(abs(ld_obj$pos - lead_pos))
  v <- as.numeric(ld_obj$mat[lead_row, ])^2; v[v > 1] <- 1
  dt <- data.table(position = as.integer(ld_obj$pos), r2 = v)
  dt[, .(r2 = max(r2, na.rm = TRUE)), by = position]
}
cat("Loading EUR LD and computing r^2 to lead...\n")
eur_ld <- load_ld_for_window(LD_BASE, LD_BLOCKS)
eur_r2 <- r2_to_lead(eur_ld, EUR_LEAD_POS)
eur_locus <- merge(eur_locus, eur_r2, by = "position", all.x = TRUE)
eur_locus[is.na(r2), r2 := 0]
cat(sprintf("  LD merged: %d / %d EUR variants matched\n",
            sum(eur_locus$r2 > 0), nrow(eur_locus)))

# ── Inline SuSiE fallback (when no stored PIPs for this study/window) ─────────
# Mirrors figS09_locus_zoom_pg.R's inline UKBB-EUR block: L=1 single-effect fit
# on the EUR z-scores + 1kg LD, conservative for these clean single-signal loci.
if (!nrow(pips) && !is.null(eur_ld) && requireNamespace("susieR", quietly = TRUE)) {
  cat("  Running inline SuSiE on EUR sumstats (no stored PIPs)...\n")
  shared_pos <- intersect(eur_locus$position, as.integer(rownames(eur_ld$mat)))
  if (length(shared_pos) >= 30) {
    sub <- eur_locus[position %in% shared_pos][order(position)]
    pc  <- as.character(sub$position); R <- eur_ld$mat[pc, pc]
    z   <- sub$beta / sub$se
    if (length(z) == nrow(R)) {
      fit <- tryCatch(susieR::susie_rss(z = z, R = R, L = 1, coverage = 0.95,
                        n = cfg$n_gwas, max_iter = 500, refine = FALSE),
                      error = function(e) { cat("  susie_rss failed:", e$message, "\n"); NULL })
      if (!is.null(fit) && (isTRUE(fit$converged) || max(fit$pip) >= 0.3)) {
        cs_idx <- integer(length(fit$pip))
        if (length(fit$sets$cs)) for (k in seq_along(fit$sets$cs)) cs_idx[fit$sets$cs[[k]]] <- k
        pips <- data.table(study = STUDY, chromosome = CHR, position = sub$position,
                           susie_pip = fit$pip, susie_cs = cs_idx, susie_converged = TRUE)[susie_pip > 0]
        pip_source <- "inline"
        cat(sprintf("  Inline SuSiE: %d variants (max PIP=%.3f)\n", nrow(pips), max(pips$susie_pip)))
      } else cat("  Inline SuSiE did not converge.\n")
    }
  }
}

# ── eQTL track ────────────────────────────────────────────────────────────────
eqtl_file  <- file.path(EQTL_DIR, paste0("chr", CHR, "_marginal_summary_results.tsv"))
eqtl_track <- NULL; best_gene <- EQTL_GENE
if (file.exists(eqtl_file)) {
  eqtl <- fread(eqtl_file, select = c("Variant","CHR","POS","Beta","SE","PVAL","GeneSymbol"))
  eqtl_loc <- eqtl[CHR == CHR & POS >= WIN_START & POS <= WIN_END]
  if (nrow(eqtl_loc[GeneSymbol == EQTL_GENE]) > 0) {
    eqtl_track <- eqtl_loc[GeneSymbol == EQTL_GENE]; eqtl_track[, p := as.numeric(PVAL)]
    cat(sprintf("  eQTL track: %s (%d variants)\n", best_gene, nrow(eqtl_track)))
  } else cat("  eQTL: focal gene not in window — eQTL track skipped.\n")
}

# ── SuSiE-COLOC top variant (gold diamond) ───────────────────────────────────
# coloc_pp4 = locus-level PP.H4.susie (whole-signal sharing); coloc_top_pp =
# SNP.PP.H4 of the single shared variant (top_snp_PP) — how confidently coloc
# resolves WHICH variant is shared. The latter is plotted on the PIP track: it
# is typically far sharper than the GWAS-alone SuSiE PIP (which LD smears),
# because conditioning on the eQTL concentrates the shared-variant posterior.
coloc_top_pos <- NA_integer_; coloc_pp4 <- NA_real_; coloc_abf_pp4 <- NA_real_
coloc_top_pp <- NA_real_; coloc_method <- NA_character_
coloc_gwas_ancestry <- NA_character_; coloc_ld_panel <- NA_character_
coloc_ld_panel_n <- NA_real_; coloc_ld_reliability <- NA_character_
coloc_locus_lambda_s <- NA_real_; coloc_provenance <- NA_character_

# Main Figure 2 must read the synchronized promoted aggregate, not a per-study
# method directory that can predate promotion. The latter remains a legacy
# fallback only for general-purpose supplementary renders.
if (file.exists(COLOC_INPUT)) {
  required_coloc <- c(
    "gwas_name", "gene", "PP.H4.susie", "PP.H4.abf", "top_snp",
    "top_snp_PP", "method", "ancestry", "ld_panel", "ld_panel_n",
    "ld_reliability", "locus_lambda_s"
  )
  tt <- fread(COLOC_INPUT, select = required_coloc)
  if (!all(required_coloc %in% names(tt))) {
    stop("Promoted COLOC table lacks required columns: ",
         paste(setdiff(required_coloc, names(tt)), collapse = ", "))
  }
  hit <- tt[
    gwas_name == STUDY & gene == EQTL_GENE &
      is.finite(PP.H4.susie) & PP.H4.susie > 0.5
  ]
  if (MATCHED_MAIN && nrow(hit) != 1L) {
    stop("Expected exactly one promoted multi-signal row for ", EQTL_GENE,
         " x ", STUDY, "; found ", nrow(hit))
  }
  if (nrow(hit)) {
    setorder(hit, -PP.H4.susie, -top_snp_PP)
    hit <- hit[1]
    pos <- as.integer(sub("^[0-9]+:", "", hit$top_snp))
    if (!is.na(pos)) {
      coloc_top_pos <- pos
      coloc_pp4 <- hit$PP.H4.susie
      coloc_abf_pp4 <- hit$PP.H4.abf
      coloc_top_pp <- hit$top_snp_PP
      coloc_method <- hit$method
      coloc_gwas_ancestry <- hit$ancestry
      coloc_ld_panel <- hit$ld_panel
      coloc_ld_panel_n <- hit$ld_panel_n
      coloc_ld_reliability <- hit$ld_reliability
      coloc_locus_lambda_s <- hit$locus_lambda_s
      coloc_provenance <- normalizePath(COLOC_INPUT)
      cat(sprintf(
        "  COLOC top SNP (%s x %s, promoted aggregate): chr%d:%d  PP4=%.3f  SNP.PP.H4=%.3f\n",
        EQTL_GENE, STUDY, CHR, coloc_top_pos, coloc_pp4, coloc_top_pp
      ))
    }
  }
}
if (is.na(coloc_top_pos) && MATCHED_MAIN) {
  stop("No promoted multi-signal COLOC row for ", EQTL_GENE, " x ", STUDY,
       " in ", COLOC_INPUT)
}
if (is.na(coloc_top_pos)) {
  for (d in c("susie_coloc", "susie_coloc_polyfun", "susie_coloc_1kg")) {
    f <- file.path(BASE, "results", d, STUDY, sprintf("susie_coloc_chr%d.csv", CHR))
    if (!file.exists(f)) next
    tt <- fread(f)
    hit <- tt[gene == EQTL_GENE & method == "susie" &
                is.finite(PP.H4.susie) & PP.H4.susie > 0.5]
    if (!nrow(hit)) next
    setorder(hit, -top_snp_PP)
    pos <- as.integer(sub("^[0-9]+:", "", hit$top_snp[1]))
    if (!is.na(pos)) {
      coloc_top_pos <- pos
      coloc_pp4 <- hit$PP.H4.susie[1]
      coloc_abf_pp4 <- hit$PP.H4.abf[1]
      coloc_top_pp <- hit$top_snp_PP[1]
      coloc_method <- hit$method[1]
      coloc_provenance <- normalizePath(f)
      cat(sprintf(
        "  COLOC top SNP (%s x %s, legacy %s): chr%d:%d  PP4=%.3f  SNP.PP.H4=%.3f\n",
        EQTL_GENE, STUDY, d, CHR, coloc_top_pos, coloc_pp4, coloc_top_pp
      ))
      break
    }
  }
}

# ── Atlas logFC for gene-track coloring ──────────────────────────────────────
LFC_SAT <- 0.3; LFC_MIN <- 0.05   # saturate at the Tier-1 DEG threshold (|log2FC|>=0.3)
                                  # so real DEGs render at full color instead of near-white
gene_atlas <- if (file.exists(ATLAS_FILE))
  fread(ATLAS_FILE, select = c("human_symbol","bulk_logFC","bulk_padj")) else
  data.table(human_symbol = character(), bulk_logFC = double(), bulk_padj = double())
# Diverging: deep blue (down) -> white -> deep magenta (up) — house semantic
# (up = #C9265E, down = #1565C0; was red #C62828, 2026-07-01 fix).
deg_palette <- colorRampPalette(c("#1565C0","#90CAF9","#FFFFFF","#F48FB1","#C9265E"))(101)
lfc_to_color <- function(lfc, padj) {
  if (is.na(lfc) || is.na(padj) || padj >= 0.05 || abs(lfc) < LFC_MIN) return("#9E9E9E")
  v <- pmin(pmax(lfc / LFC_SAT, -1), 1); deg_palette[round((v + 1) * 50) + 1]
}

# ─────────────────────────────────────────────────────────────────────────────
# Render
# ─────────────────────────────────────────────────────────────────────────────
cat("Output:", out_path, "\n")
# Compact main-text geometry (2026-07-06): ~4.4 x ~4.66 in (was 5.4 x 4.21) so the
# 6pt font reads proportionally next to the other Fig 2 panels. Narrower page +
# tightened per-track heights/gaps; MARGIN_R held at 1.00 for the right-margin legends.
# Sized to the Figure 2 contract (figure2_panel_sizes.tsv), so Fig2G/H
# place at 100% next to the other Fig2 panels. Tracks compressed ~0.69x vs the 4.4x4.66 version.
MARGIN_L <- if (MATCHED_MAIN) 0.68 else 0.62
# Right margin holds the shared legend and/or the track-type labels. Shrink it when
# those are stripped (F: neither; G: labels only) so the plot expands into the freed
# width and the page tightens naturally. Plot region (PLOT_W) is constant so F/G tracks align.
MARGIN_R <- if (MATCHED_MAIN) 0.10 else if (NO_LEGEND && NO_TRACKLABEL) 0.15 else if (NO_LEGEND) 0.55 else 1.00
PLOT_W  <- if (MATCHED_MAIN) 2.48 else 2.33
PLOT_X  <- MARGIN_L; PAGE_W <- MARGIN_L + PLOT_W + MARGIN_R
# track layout (inches from top); font floor 6pt
if (MATCHED_MAIN) {
  # Exact 3.26 x 3.10-in main-panel geometry. Height is removed from the PIP
  # track and inter-track gaps; the GWAS and eQTL clouds retain their heights.
  y1 <- 0.10; h1 <- 0.76
  y3 <- y1 + h1 + 0.08; h3 <- 0.38
  y4 <- y3 + h3 + 0.06; h4 <- 0.62
  y5 <- y4 + h4 + 0.08; h5 <- 0.54
  y6 <- y5 + h5 + 0.04
  PAGE_H <- 3.10
} else {
  y1 <- 0.32; h1 <- 0.76                       # EUR GWAS (top room for the per-panel header)
  y3 <- y1 + h1 + 0.12; h3 <- 0.62             # PIP
  y4 <- y3 + h3 + 0.09; h4 <- 0.62             # eQTL (taller -> scatter spreads out)
  y5 <- y4 + h4 + 0.12; h5 <- 0.54             # genes (room for 2-tier de-collided labels)
  y6 <- y5 + h5 + 0.05                         # ruler
  PAGE_H <- y6 + 0.30
}

cairo_pdf(out_path, width = PAGE_W, height = PAGE_H, family = "Helvetica")   # embeds real Helvetica (plain pdf() did not; family= needed so grid text isn't Nimbus)
pageCreate(width = PAGE_W, height = PAGE_H, default.units = "inches",
           showGuides = FALSE, xgrid = 0, ygrid = 0)
top_to_grid <- function(y_top) PAGE_H - y_top

# The compact main panels are title-free for publication. Locus identity remains
# in the caption and source sidecar. General-purpose outputs retain their
# historical inline header.
if (!MATCHED_MAIN) local({
  hg <- EQTL_GENE
  hr <- paste0("  ·  ", TRAIT_TOKEN, " (EUR)")
  plotText(hg, x = PLOT_X, y = 0.15, fontsize = 6, fontface = "italic",
           just = c("left", "center"), default.units = "inches")
  gw <- convertWidth(grobWidth(textGrob(hg, gp = gpar(fontsize = 6, fontface = "italic"))),
                     "inches", valueOnly = TRUE)
  plotText(hr, x = PLOT_X + gw, y = 0.15, fontsize = 6,
           just = c("left", "center"), default.units = "inches")
})

# Single-hue sequential gray->navy gradient (2026-07-01, replaces the rainbow
# grey/blue/green/yellow/orange/red ramp per FIGURE_GUIDELINES.md "Never
# rainbow/jet" — matches figS09_locus_zoom_pg.R's enzyme-panel LD scale).
ld_palette <- colorRampPalette(c("#E6E6E5", "#1B2F5B"))
ld_pal_vec <- ld_palette(100)

track_label <- function(label, x = PAGE_W - MARGIN_R + 0.05, y)
  plotText(label = label, x = x, y = y, fontsize = 6,
           just = c("left","top"), default.units = "inches")
axis_label <- function(label, y, h) {
  if (NO_YTITLE) return(invisible())   # right panel shares the left panel's y-axis titles
  xoff <- if (MATCHED_MAIN) 0.49 else 0.42
  plotText(label = label, rot = 90, x = MARGIN_L - xoff, y = y + h/2,
           fontsize = 6, just = "center", default.units = "inches")
}
# Draws one point + label legend item at (x, y) inches-from-top and returns the
# x position the NEXT item should start at, measuring the label's actual
# rendered width (grid::grobWidth) instead of a hand-guessed offset — this is
# what a fixed offset got wrong (2026-07-01: legend items were overlapping).
legend_item <- function(x, y, col, label, fill = col, pch = 19, size = 0.07, gap = 0.14) {
  grid.points(unit(x, "inches"), unit(top_to_grid(y), "inches"),
              pch = pch, size = unit(size, "inches"), gp = gpar(col = col, fill = fill, lwd = 0.5))
  lab_x <- x + size / 2 + 0.05
  plotText(label, x = lab_x, y = y, fontsize = 6, fontcolor = "grey25",
           just = c("left", "center"), default.units = "inches")
  txt_w <- convertWidth(grobWidth(textGrob(label, gp = gpar(fontsize = 6))),
                        "inches", valueOnly = TRUE)
  lab_x + txt_w + gap
}
draw_axis_spines <- function(y, h) {
  yt <- top_to_grid(y); yb <- top_to_grid(y + h)
  grid.segments(unit(PLOT_X,"inches"), unit(yb,"inches"), unit(PLOT_X,"inches"), unit(yt,"inches"),
                gp = gpar(col = "black", lwd = 0.5))
  grid.segments(unit(PLOT_X,"inches"), unit(yb,"inches"), unit(PLOT_X+PLOT_W,"inches"), unit(yb,"inches"),
                gp = gpar(col = "black", lwd = 0.5))
}

message(sprintf("[caption] %s  ·  %s", GENE, cfg$gwas_label))

# Track 1: EUR GWAS, LD-colored ───────────────────────────────────────────────
eur_pg <- eur_locus[!is.na(p), .(chrom = paste0("chr", CHR), pos = position, p = p, r2 = r2)]
setorder(eur_pg, r2)
eur_pg[, fill_col := ld_pal_vec[pmin(100, pmax(1, ceiling(r2 * 100)))]]
eur_pg_df <- as.data.frame(eur_pg)
eur_y_max <- ceiling(max(-log10(eur_pg_df$p), na.rm = TRUE)) + 1
mh_eur <- plotManhattan(data = eur_pg_df, chrom = paste0("chr", CHR),
  chromstart = WIN_START, chromend = WIN_END, assembly = "hg19",
  fill = NA, pch = 19, cex = 0.001, sigVal = 5e-8, sigLine = TRUE, sigCol = "black",
  range = c(0, eur_y_max), x = PLOT_X, width = PLOT_W, y = y1, height = h1,
  default.units = "inches")
pos_to_x <- function(pos) PLOT_X + (pos - WIN_START) / (WIN_END - WIN_START) * PLOT_W
eur_p_to_y <- function(pv) { top <- PAGE_H - y1; bot <- PAGE_H - (y1 + h1)
  bot + (-log10(pv)) / eur_y_max * (top - bot) }
grid.points(unit(pos_to_x(eur_pg_df$pos),"inches"), unit(eur_p_to_y(eur_pg_df$p),"inches"),
            pch = 19, size = unit(0.03,"inches"), gp = gpar(col = eur_pg_df$fill_col))
# lead diamond (signal guaranteed: these are colocalizing loci)
lead_p <- eur_pg_df$p[which.min(abs(eur_pg_df$pos - EUR_LEAD_POS))]
grid.points(unit(pos_to_x(EUR_LEAD_POS),"inches"), unit(eur_p_to_y(lead_p),"inches"),
            pch = 23, size = unit(0.085,"inches"), gp = gpar(col = "black", fill = "#C2185B", lwd = 0.8))
annoYaxis(plot = mh_eur, at = pretty(c(0, max(-log10(eur_pg$p)))), fontsize = 6)
draw_axis_spines(y1, h1); axis_label(if (MATCHED_MAIN) "GWAS −log10 p" else "-log10(p)", y1, h1)
if (!NO_TRACKLABEL) track_label("GWAS", y = y1 + 0.02)   # generic track type; locus/trait is in the header. Shared once on left panel.

# shared LD legend (right margin) — drawn once across the D|E pair (skipped on the left panel)
if (!NO_LEGEND) {
ld_x0 <- PAGE_W - MARGIN_R + 0.05; ld_y <- y1 + 0.32; ld_w <- 0.85
for (k in 0:49) grid.rect(unit(ld_x0 + k*(ld_w/50),"inches"), unit(top_to_grid(ld_y),"inches"),
  width = unit(ld_w/50 + 0.005,"inches"), height = unit(0.10,"inches"), just = c("left","center"),
  gp = gpar(col = NA, fill = ld_pal_vec[pmin(100, pmax(1, ceiling((k/49)*100)))]))
plotText("0.0", x = ld_x0, y = ld_y + 0.10, fontsize = 6, just = c("left","top"), default.units = "inches")
plotText("0.5", x = ld_x0 + ld_w/2, y = ld_y + 0.10, fontsize = 6, just = c("center","top"), default.units = "inches")
plotText("1.0", x = ld_x0 + ld_w, y = ld_y + 0.10, fontsize = 6, just = c("right","top"), default.units = "inches")
plotText("r² to lead", x = ld_x0 + ld_w/2, y = ld_y - 0.10, fontsize = 6,
         just = c("center","bottom"), default.units = "inches")
grid.points(unit(ld_x0 + 0.05,"inches"), unit(top_to_grid(ld_y + 0.30),"inches"),
            pch = 23, size = unit(0.10,"inches"), gp = gpar(col = "black", fill = "#C2185B", lwd = 0.6))
plotText("Lead SNP", x = ld_x0 + 0.15, y = ld_y + 0.30, fontsize = 6, just = c("left","center"), default.units = "inches")
}

# Track 2: SuSiE PIP lollipops (EUR) ──────────────────────────────────────────
pip_top_t <- y3 + 0.10; pip_base_t <- y3 + h3 - 0.10
pip_base <- top_to_grid(pip_base_t); pip_top <- top_to_grid(pip_top_t)
pip_to_y <- function(p) pip_base + p * (pip_top - pip_base)
grid.segments(unit(PLOT_X,"inches"), unit(pip_base,"inches"), unit(PLOT_X+PLOT_W,"inches"), unit(pip_base,"inches"),
              gp = gpar(col = "black", lwd = 0.6))
grid.segments(unit(PLOT_X,"inches"), unit(pip_base,"inches"), unit(PLOT_X,"inches"), unit(pip_top,"inches"),
              gp = gpar(col = "black", lwd = 0.6))
grid.segments(unit(PLOT_X,"inches"), unit(pip_to_y(0.9),"inches"), unit(PLOT_X+PLOT_W,"inches"), unit(pip_to_y(0.9),"inches"),
              gp = gpar(col = "grey50", lty = 2, lwd = 0.4))
for (val in c(0, 0.5, 1.0)) {
  grid.segments(unit(PLOT_X-0.05,"inches"), unit(pip_to_y(val),"inches"), unit(PLOT_X,"inches"), unit(pip_to_y(val),"inches"),
                gp = gpar(col = "black", lwd = 0.5))
  grid.text(format(val, nsmall = 1), unit(PLOT_X-0.08,"inches"), unit(pip_to_y(val),"inches"),
            just = "right", gp = gpar(fontsize = 6))
}
if (nrow(pips)) for (i in seq_len(nrow(pips))) {
  xi <- pos_to_x(pips$position[i]); yi <- pip_to_y(pips$susie_pip[i])
  in_cs <- !is.na(pips$susie_cs[i]) && pips$susie_cs[i] > 0
  grid.segments(unit(xi,"inches"), unit(pip_base,"inches"), unit(xi,"inches"), unit(yi,"inches"),
                gp = gpar(col = "grey85", lwd = 0.3))
  grid.points(unit(xi,"inches"), unit(yi,"inches"), pch = if (in_cs) 19 else 21,
              size = unit(if (in_cs) 0.055 else 0.028,"inches"),
              gp = gpar(col = "#1565C0", fill = if (in_cs) "#1565C0" else "white", lwd = 0.5))
}
# Gold lollipop = coloc-resolved shared variant (SNP.PP.H4). Sharper than the
# GWAS-alone PIP because it conditions on the (well-resolved) eQTL signal.
if (!is.na(coloc_top_pos) && !is.na(coloc_top_pp) &&
    coloc_top_pos >= WIN_START && coloc_top_pos <= WIN_END) {
  cx <- pos_to_x(coloc_top_pos); cy <- pip_to_y(coloc_top_pp)
  grid.segments(unit(cx,"inches"), unit(pip_base,"inches"), unit(cx,"inches"), unit(cy,"inches"),
                gp = gpar(col = "#F9A825", lwd = 0.9))
  grid.points(unit(cx,"inches"), unit(cy,"inches"), pch = 23, size = unit(0.085,"inches"),
              gp = gpar(col = "black", fill = "#FFD600", lwd = 0.8))
  lab_left <- coloc_top_pos > (WIN_START + 0.6 * (WIN_END - WIN_START))
  plotText(sprintf("coloc shared variant\nSNP.PP.H4 = %.2f", coloc_top_pp),
           x = if (MATCHED_MAIN) PLOT_X + PLOT_W - 0.03 else if (lab_left) cx - 0.10 else cx + 0.10,
           y = if (MATCHED_MAIN) y3 + 0.03 else y3 + 0.14,
           fontsize = 6, just = c(if (MATCHED_MAIN || lab_left) "right" else "left","top"),
           default.units = "inches", fontcolor = "#8D6E00")
}
axis_label(if (MATCHED_MAIN) "SuSiE PIP" else "PIP", y3, h3)
if (!NO_TRACKLABEL) track_label("SuSiE PIP", y = y3 + 0.02)
# mini legend — width-measured (legend_item()), not hand-guessed offsets. Shared across D|E (skip on left panel).
if (!NO_LEGEND) {
xc <- legend_item(PLOT_X + 0.10, y3 + 0.10, "#1565C0", "EUR (in CS)")
xc <- legend_item(xc, y3 + 0.10, "#9E9E9E", "not in CS", fill = "white", pch = 21, size = 0.05)
xc <- legend_item(xc, y3 + 0.10, "black", "coloc shared variant", fill = "#FFD600", pch = 23, size = 0.08)
}

# Track 3: eQTL Manhattan ─────────────────────────────────────────────────────
if (!is.null(eqtl_track) && nrow(eqtl_track) > 0) {
  eqtl_pg <- eqtl_track[!is.na(p) & p > 0, .(pos = POS, p = p)]
  eqtl_y_max <- ceiling(max(-log10(eqtl_pg$p), na.rm = TRUE)) + 1
  mh_eqtl <- plotManhattan(data = data.frame(chrom = paste0("chr", CHR), pos = eqtl_pg$pos, p = eqtl_pg$p),
    chrom = paste0("chr", CHR), chromstart = WIN_START, chromend = WIN_END, assembly = "hg19",
    fill = "#E07B39", pch = 19, cex = 0.16, sigVal = 1e-5, sigLine = FALSE,
    range = c(0, eqtl_y_max), x = PLOT_X, width = PLOT_W, y = y4, height = h4, default.units = "inches")
  annoYaxis(plot = mh_eqtl, at = pretty(c(0, max(-log10(eqtl_pg$p)))), fontsize = 6)
  draw_axis_spines(y4, h4); axis_label(if (MATCHED_MAIN) "eQTL −log10 p" else "-log10(p)", y4, h4)
  if (!NO_TRACKLABEL) track_label("eQTL", y = y4 + 0.02)   # generic; focal gene is in the header + pink in the gene track
  if (!is.na(coloc_top_pos) && coloc_top_pos >= WIN_START && coloc_top_pos <= WIN_END) {
    coloc_p <- eqtl_pg$p[which.min(abs(eqtl_pg$pos - coloc_top_pos))]
    cyt <- PAGE_H - y4; cyb <- PAGE_H - (y4 + h4)
    cy <- cyb + (-log10(coloc_p)) / eqtl_y_max * (cyt - cyb)
    grid.points(unit(pos_to_x(coloc_top_pos),"inches"), unit(cy,"inches"),
                pch = 23, size = unit(0.085,"inches"), gp = gpar(col = "black", fill = "#FFD600", lwd = 0.8))
    plotText(sprintf("COLOC top SNP\nPP4 = %.2f", coloc_pp4),
             x = if (MATCHED_MAIN) PLOT_X + PLOT_W - 0.03 else pos_to_x(coloc_top_pos) + 0.10,
             y = if (MATCHED_MAIN) y4 + 0.03 else y4 + 0.10,
             fontsize = 6, just = c(if (MATCHED_MAIN) "right" else "left","top"),
             default.units = "inches", fontcolor = "grey20")
  }
} else {
  plotText(sprintf("(no %s eQTL in window)", EQTL_GENE), x = PLOT_X + PLOT_W/2, y = y4 + h4/2,
           fontsize = 6, fontcolor = "grey50", default.units = "inches")
  if (!NO_TRACKLABEL) track_label("eQTL", y = y4 + 0.02)   # generic; focal gene is in the header + pink in the gene track
}

# Track 4: gene track (custom transcript model) ──────────────────────────────
locus_gr <- GRanges(paste0("chr", CHR), IRanges(WIN_START, WIN_END))
genes_in_win <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                  filter = list(tx_chrom = paste0("chr", CHR))))
genes_in_win <- subsetByOverlaps(genes_in_win, locus_gr)
gene_symbols <- if (length(genes_in_win) > 0) unique(na.omit(suppressMessages(
  mapIds(org.Hs.eg.db, keys = genes_in_win$gene_id, column = "SYMBOL",
         keytype = "ENTREZID", multiVals = "first")))) else character(0)
gene_colors <- setNames(character(length(gene_symbols)), gene_symbols)
for (g in gene_symbols) { hit <- gene_atlas[human_symbol == g][1]
  gene_colors[g] <- if (nrow(hit) && !is.na(hit$bulk_logFC)) lfc_to_color(hit$bulk_logFC, hit$bulk_padj) else "#9E9E9E" }
xin <- function(p) PLOT_X + (pmax(WIN_START, pmin(WIN_END, p)) - WIN_START) / (WIN_END - WIN_START) * PLOT_W
ex_by_gene <- suppressMessages(exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene"))
gmod <- list()
for (i in seq_along(genes_in_win)) {
  gid <- as.character(genes_in_win$gene_id[i])
  sym <- suppressMessages(mapIds(org.Hs.eg.db, gid, "SYMBOL", "ENTREZID"))
  if (is.na(sym) || is.null(ex_by_gene[[gid]])) next
  exr <- reduce(ex_by_gene[[gid]])
  if (max(end(exr)) < WIN_START || min(start(exr)) > WIN_END) next
  gmod[[length(gmod)+1]] <- list(sym = sym, strand = as.character(strand(exr))[1],
    gs = min(start(exr)), ge = max(end(exr)), es = start(exr), ee = end(exr),
    col = if (sym %in% names(gene_colors)) gene_colors[[sym]] else "#9E9E9E")
}
if (length(gmod)) {
  gmod <- gmod[order(sapply(gmod, `[[`, "gs"))]
  row_end <- numeric(0)
  for (j in seq_along(gmod)) {
    gs <- max(WIN_START, gmod[[j]]$gs); placed <- FALSE
    for (r in seq_along(row_end)) if (gs > row_end[r] + 0.02*(WIN_END-WIN_START)) {
      gmod[[j]]$row <- r; row_end[r] <- min(WIN_END, gmod[[j]]$ge); placed <- TRUE; break }
    if (!placed) { row_end <- c(row_end, min(WIN_END, gmod[[j]]$ge)); gmod[[j]]$row <- length(row_end) }
  }
  nrows <- max(length(row_end), 1); row_gap <- min(0.30, (h5 - 0.20) / nrows)
  glab <- list()
  for (g in gmod) {
    yc <- y5 + 0.18 + (g$row - 1) * row_gap; yg <- top_to_grid(yc)
    x0 <- xin(g$gs); x1 <- xin(g$ge)
    grid.lines(unit(c(x0, x1),"inches"), unit(yg,"inches"), gp = gpar(col = g$col, lwd = 1.0))
    nch <- floor((x1 - x0) / 0.22)
    if (nch >= 1) { chx <- seq(x0 + 0.06, x1 - 0.06, length.out = nch + 1)
      for (cx in chx) grid.text(if (g$strand == "-") "<" else ">", unit(cx,"inches"), unit(yg,"inches"),
        gp = gpar(col = g$col, fontsize = 6)) }
    for (k in seq_along(g$es)) { ex0 <- xin(g$es[k]); ex1 <- xin(g$ee[k])
      grid.rect(unit(ex0,"inches"), unit(yg,"inches"), width = unit(max(ex1-ex0, 0.006),"inches"),
        height = unit(0.085,"inches"), just = c("left","center"), gp = gpar(col = NA, fill = g$col)) }
    glab[[length(glab)+1]] <- list(sym = g$sym, xmid = (x0+x1)/2, col = g$col, yc = yc, row = g$row)
  }
  # De-collide gene-symbol labels in three global vertical lanes. Per-body-row
  # tiers can coincide across adjacent rows (the former MIR4780/FABP1 collision),
  # so lane occupancy must be tracked across the complete gene track.
  for (i in seq_along(glab)) glab[[i]]$w <- convertWidth(grobWidth(textGrob(
    glab[[i]]$sym, gp = gpar(fontsize = 6, fontface = "italic"))), "inches", valueOnly = TRUE)
  glab <- glab[order(sapply(glab, `[[`, "xmid"))]
  LABEL_Y <- y5 + c(0.06, 0.265, 0.48)
  redge <- rep(-Inf, length(LABEL_Y))
  for (L in glab) {
    lft <- L$xmid - L$w/2
    available <- which(lft > redge + 0.02)
    t <- if (length(available)) available[1] else which.min(redge)
    redge[t] <- L$xmid + L$w/2
    grid.text(L$sym, unit(L$xmid,"inches"), unit(top_to_grid(LABEL_Y[t]),"inches"),
              gp = gpar(col = L$col, fontsize = 6, fontface = "italic"))
  }
}
if (MATCHED_MAIN) axis_label("Genes", y5, h5)
if (!NO_TRACKLABEL) track_label("Genes", y = y5 + 0.02)
# logFC scale bar — shared legend, drawn once across the D|E pair (skipped on the left panel)
if (!NO_LEGEND) {
sb_x0 <- PAGE_W - MARGIN_R + 0.05; sb_y <- y5 + 0.32; sb_w <- 0.85
for (k in 0:49) { step_lfc <- (-LFC_SAT) + (2*LFC_SAT)*(k/49)
  grid.rect(unit(sb_x0 + k*(sb_w/50),"inches"), unit(top_to_grid(sb_y),"inches"),
    width = unit(sb_w/50 + 0.005,"inches"), height = unit(0.10,"inches"), just = c("left","center"),
    gp = gpar(col = NA, fill = deg_palette[round((step_lfc/LFC_SAT + 1)*50) + 1])) }
plotText(sprintf("-%.1f", LFC_SAT), x = sb_x0, y = sb_y + 0.10, fontsize = 6, just = c("left","top"), default.units = "inches")
plotText("0", x = sb_x0 + sb_w/2, y = sb_y + 0.10, fontsize = 6, just = c("center","top"), default.units = "inches")
plotText(sprintf("+%.1f", LFC_SAT), x = sb_x0 + sb_w, y = sb_y + 0.10, fontsize = 6, just = c("right","top"), default.units = "inches")
plotText("log2FC", x = sb_x0 + sb_w/2, y = sb_y - 0.10, fontsize = 6, just = c("center","bottom"), default.units = "inches")
grid.rect(unit(sb_x0,"inches"), unit(top_to_grid(sb_y + 0.30),"inches"), width = unit(0.10,"inches"),
          height = unit(0.10,"inches"), just = c("left","center"), gp = gpar(col = NA, fill = "#9E9E9E"))
plotText("n.s.", x = sb_x0 + 0.13, y = sb_y + 0.30, fontsize = 6, just = c("left","center"), default.units = "inches")
}

# Track 5: genome ruler ───────────────────────────────────────────────────────
plotGenomeLabel(chrom = paste0("chr", CHR), chromstart = WIN_START, chromend = WIN_END,
  assembly = "hg19", scale = "Mb", commas = TRUE, sequence = FALSE, fontsize = 6,
  x = PLOT_X, y = y6, length = PLOT_W, default.units = "inches")

# Guides are already disabled in pageCreate(). Calling pageGuideHide() here
# invokes a grid removal pass that plotgardener documents as creating a second
# PDF page, so close the device directly.
dev.off()
normalize_page(out_path, PAGE_W, PAGE_H)

# ── sidecar source CSV ────────────────────────────────────────────────────────
src <- data.table(gene = GENE, eqtl_gene = EQTL_GENE, gwas = STUDY, trait = cfg$trait,
  chr = CHR, lead_pos_hg19 = LEAD_POS, win_start = WIN_START, win_end = WIN_END,
  coloc_top_snp_pos = coloc_top_pos, coloc_PP4_susie = round(coloc_pp4, 4),
  coloc_PP4_abf = round(coloc_abf_pp4, 4), coloc_method = coloc_method,
  coloc_shared_variant_PP = round(coloc_top_pp, 4),
  coloc_input = coloc_provenance, gwas_ancestry = coloc_gwas_ancestry,
  coloc_ld_panel = coloc_ld_panel, coloc_ld_panel_n = coloc_ld_panel_n,
  coloc_ld_reliability = coloc_ld_reliability,
  coloc_locus_lambda_s = coloc_locus_lambda_s,
  eqtl_source = "Broadaway liver eQTL", eqtl_ancestry = "EUR",
  plot_ld_source = normalizePath(LD_BASE),
  n_pip_variants = nrow(pips), pip_source = pip_source)
fwrite(src, sub("\\.pdf$", "_source.csv", out_path))
cat("\nDone:", out_path, "\n")
