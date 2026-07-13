#!/usr/bin/env Rscript
# ============================================================================
# q2_cutoff_bottleneck_grid.R   (Vacca 2024 benchmark follow-up — Q2 KEY)
#
# QUESTION: the prior conserved-core (1,355 genes) ∩ Vacca-951 overlap (172,
# Fisher OR 2.8x) is bottlenecked by ARBITRARY human + mouse LFC/padj cutoffs.
# Quantify how (i) core SIZE, (ii) core∩Vacca-951 OVERLAP, (iii) Fisher OR move
# across a GRID of selection thresholds, to show 1355/172/2.8x is ONE arbitrary
# point on a continuum.
#
# We RECONSTRUCT the per-gene concordance from raw inputs, mirroring the canonical
# builder 01_corrected_gene_concordance.R EXACTLY (nafl_vs_nash human anchor, 4
# mouse diets via orthologs, "concordant in diet" = h_padj<P_H & m_padj<P_M &
# sign(h_lfc)==sign(m_lfc), "Conserved core" = n_concordant >= N_REQ), but we
# PARAMETERIZE the four cutoffs and additionally ADD a human |logFC| floor (the
# canonical pipeline applies NO LFC filter to the concordance call — only padj +
# sign — so |logFC|=0 reproduces canonical).
#
# GRID:
#   human |logFC|  in {0.25, 0.5, 0.75, 1.0}   (+ 0.0 = canonical, no LFC filter)
#   human padj     in {0.01, 0.05, 0.1}
#   mouse padj     in {0.05, 0.1}
#   n_concordant   in {2, 3, 4}
#
# Vacca-951 = MOESM4 Table S4, Gene_used_in_DSEA==1. Fair universe = our paired
# concordance genes ∩ Vacca S4 tested universe (same as script 07).
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(readxl) })

BASE  <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
VACCA <- file.path(BASE, "data/external/vacca_2024")
H_INT <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
ANNOT <- file.path(H_INT, "results/gene_annotation")
DS_DIR<- file.path(H_INT, "results/disease_signatures")
MPD   <- file.path(BASE, "RNA-seq/Mouse/Unified_Integration/results/per_diet")
OUT   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
strip <- function(x) gsub("\\..*", "", x)

DIETS <- c("MCD", "HFD", "CDAHFD", "FPC")

cat(strrep("=", 78), "\n")
cat("Q2 cutoff-bottleneck grid: conserved-core size / Vacca-951 overlap / Fisher OR\n")
cat(strrep("=", 78), "\n")

# ── Ortholog map ─────────────────────────────────────────────────────────────
ortho <- fread(file.path(ANNOT, "ortholog_mapping.tsv"))

# ── Human nafl_vs_nash anchor (canonical builder uses nafl_vs_nash_dream.csv;
#    we mirror that EXACTLY for canonical reproduction) ────────────────────────
hraw <- fread(file.path(DS_DIR, "nafl_vs_nash_dream.csv"))
hraw[, gene_base := strip(gene)]
# padj column
if (!"padj" %in% names(hraw) && "adj.P.Val" %in% names(hraw)) setnames(hraw, "adj.P.Val", "padj")
h_mapped <- merge(
  hraw[, .(gene_base, h_lfc = logFC, h_padj = padj)],
  ortho[, .(human_gene_id, mouse_gene_id, human_symbol)],
  by.x = "gene_base", by.y = "human_gene_id")

# ── 4 mouse diets, long-format pairing (one row per gene×diet) ────────────────
mouse_long <- rbindlist(lapply(DIETS, function(d) {
  md <- fread(file.path(MPD, paste0(d, "_de_results.csv")))
  md[, mouse_base := strip(gene)]
  md[, .(mouse_gene_id = mouse_base, m_lfc = logFC, m_padj = adj.P.Val, diet = d)]
}))

paired <- merge(h_mapped, mouse_long, by = "mouse_gene_id", allow.cartesian = TRUE)
cat(sprintf("Paired human-anchor × mouse-diet rows: %d (genes: %d)\n",
            nrow(paired), uniqueN(paired$human_symbol)))

# ── Vacca-951 signature + tested universe (MOESM4 Table S4) ───────────────────
s4 <- as.data.table(read_excel(file.path(VACCA, "42255_2024_1043_MOESM4_ESM.xlsx"),
                               sheet = "Table S4"))
setnames(s4, c("GeneSymbol", "Gene_used_in_DSEA(0:No/1:Yes)"),
         c("gene", "in_dsea"), skip_absent = TRUE)
s4 <- s4[!is.na(gene) & gene != ""][!duplicated(gene)]
vacca_sig <- s4[in_dsea == 1]$gene
cat(sprintf("Vacca-951 signature: %d genes ; Vacca S4 tested universe: %d\n",
            length(vacca_sig), nrow(s4)))

# Fair universe = our paired genes ∩ Vacca S4 universe (matches script 07)
our_genes <- unique(paired$human_symbol)
universe  <- intersect(our_genes, s4$gene)
vacca_sig_u <- intersect(vacca_sig, universe)
cat(sprintf("Shared universe (our paired ∩ Vacca S4): %d ; Vacca-951 in universe: %d\n\n",
            length(universe), length(vacca_sig_u)))

# ── Fisher OR helper (within the fixed universe) ─────────────────────────────
fisher_or <- function(drawn, target_u, universe) {
  drawn   <- intersect(drawn, universe)
  a  <- length(intersect(drawn, target_u)); b <- length(drawn) - a
  cc <- length(target_u) - a;               d <- length(universe) - a - b - cc
  m  <- matrix(c(a, b, cc, d), nrow = 2)
  if (any(c(a, b, cc, d) == 0)) m <- m + 0.5
  ft <- tryCatch(fisher.test(m), error = function(e) NULL)
  if (is.null(ft)) return(list(OR = NA_real_, p = NA_real_, a = a))
  list(OR = unname(ft$estimate), p = ft$p.value, a = a)
}

# ── Build conserved core for a given cutoff set ──────────────────────────────
# concordant-in-diet = h_padj<PH & m_padj<PM & |h_lfc|>=LFC & sign(h_lfc)==sign(m_lfc)
build_core <- function(LFC, PH, PM, NREQ) {
  p <- paired[h_padj < PH & m_padj < PM & abs(h_lfc) >= LFC & sign(h_lfc) == sign(m_lfc)]
  ncon <- p[, .(n_concordant = uniqueN(diet)), by = human_symbol]
  ncon[n_concordant >= NREQ, human_symbol]
}

# ── GRID ─────────────────────────────────────────────────────────────────────
grid <- CJ(LFC  = c(0.0, 0.25, 0.5, 0.75, 1.0),   # 0.0 = canonical (no LFC filter)
           PH   = c(0.01, 0.05, 0.1),
           PM   = c(0.05, 0.1),
           NREQ = c(2L, 3L, 4L))

res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i]
  core <- build_core(g$LFC, g$PH, g$PM, g$NREQ)
  core_u <- intersect(core, universe)
  ov <- length(intersect(core_u, vacca_sig_u))
  fo <- fisher_or(core_u, vacca_sig_u, universe)
  data.table(
    human_abs_logFC = g$LFC, human_padj = g$PH, mouse_padj = g$PM, n_concordant_req = g$NREQ,
    core_size       = length(core),
    core_size_in_universe = length(core_u),
    overlap_vacca951 = ov,
    pct_core_in_vacca = round(100 * ov / max(length(core_u), 1), 1),
    fisher_OR       = round(fo$OR, 3),
    fisher_p        = signif(fo$p, 3))
}))
setorder(res, human_abs_logFC, human_padj, mouse_padj, n_concordant_req)

# ── Canonical anchor (LFC=0, PH=0.05, PM=0.05, NREQ=3) — reproduce 1355/172/2.8 ─
canon <- res[human_abs_logFC == 0.0 & human_padj == 0.05 & mouse_padj == 0.05 & n_concordant_req == 3]
cat("── CANONICAL anchor (|logFC|=0, h_padj=0.05, m_padj=0.05, n_concordant>=3) ──\n")
print(canon)
cat(sprintf("  Published reference: core=1355, overlap=172, OR=2.8x.\n"))
cat(sprintf("  Reconstructed here:  core=%d, overlap=%d, OR=%.2fx  (universe=%d).\n\n",
            canon$core_size, canon$overlap_vacca951, canon$fisher_OR, length(universe)))

fwrite(res, file.path(OUT, "q2_cutoff_bottleneck_grid.csv"))
cat(sprintf("Saved full %d-cell grid → q2_cutoff_bottleneck_grid.csv\n\n", nrow(res)))

# ── Summaries: how each axis moves core-size / overlap / OR (marginals) ──────
cat("── RANGE across the full grid (excluding LFC=0 sentinel for the LFC axis stats) ──\n")
g4 <- res[human_abs_logFC > 0]   # the 4 requested LFC levels {0.25,0.5,0.75,1.0}
cat(sprintf("  core_size:       %d  →  %d   (%.0fx span)\n",
            min(res$core_size), max(res$core_size), max(res$core_size)/max(min(res$core_size),1)))
cat(sprintf("  overlap_vacca951: %d  →  %d\n", min(res$overlap_vacca951), max(res$overlap_vacca951)))
cat(sprintf("  fisher_OR:       %.2f  →  %.2f\n", min(res$fisher_OR, na.rm=TRUE), max(res$fisher_OR, na.rm=TRUE)))
cat(sprintf("  pct of core in Vacca-951: %.1f%%  →  %.1f%%\n\n",
            min(res$pct_core_in_vacca), max(res$pct_core_in_vacca)))

marg <- function(by) {
  res[, .(n_cells = .N,
          core_min = min(core_size), core_med = as.integer(median(core_size)), core_max = max(core_size),
          ov_min = min(overlap_vacca951), ov_max = max(overlap_vacca951),
          OR_min = round(min(fisher_OR, na.rm=TRUE),2), OR_med = round(median(fisher_OR, na.rm=TRUE),2),
          OR_max = round(max(fisher_OR, na.rm=TRUE),2)), by = by][order(get(by))]
}
cat("── Marginal: by human |logFC| floor ──\n");          print(marg("human_abs_logFC"))
cat("\n── Marginal: by human padj ──\n");                  print(marg("human_padj"))
cat("\n── Marginal: by mouse padj ──\n");                  print(marg("mouse_padj"))
cat("\n── Marginal: by n_concordant required ──\n");       print(marg("n_concordant_req"))

# ── Key contrast points for the report ───────────────────────────────────────
cat("\n── Selected continuum points (OR vs size/overlap tradeoff) ──\n")
pick <- res[(human_abs_logFC %in% c(0.0, 0.5, 1.0)) & human_padj %in% c(0.01, 0.05) &
            mouse_padj == 0.05 & n_concordant_req %in% c(2, 3, 4)]
print(pick[order(-fisher_OR),
           .(human_abs_logFC, human_padj, mouse_padj, n_concordant_req,
             core_size, overlap_vacca951, pct_core_in_vacca, fisher_OR, fisher_p)])

# Max-OR cell and the "loosest" cell, for the headline
best_or <- res[which.max(fisher_OR)]
loose   <- res[human_abs_logFC == 0.0 & human_padj == 0.1 & mouse_padj == 0.1 & n_concordant_req == 2]
cat(sprintf("\n  MAX Fisher OR cell: OR=%.2f at |logFC|>=%.2f, h_padj<%.2f, m_padj<%.2f, n_con>=%d (core=%d, ov=%d)\n",
            best_or$fisher_OR, best_or$human_abs_logFC, best_or$human_padj, best_or$mouse_padj,
            best_or$n_concordant_req, best_or$core_size, best_or$overlap_vacca951))
cat(sprintf("  LOOSEST cell:       OR=%.2f, core=%d, overlap=%d (|logFC|>=0, h_padj<0.1, m_padj<0.1, n_con>=2)\n",
            loose$fisher_OR, loose$core_size, loose$overlap_vacca951))

cat("\n", strrep("=", 78), "\n")
cat("DONE. The (1355 / 172 / 2.8x) point is ONE cell of a", nrow(res), "-cell continuum.\n")
cat(strrep("=", 78), "\n")
