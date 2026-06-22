#!/usr/bin/env Rscript
# coloc_deg_convergence.R  (Team B, Fig 4 validation)
# -----------------------------------------------------------------------------
# MESSAGE: genetic-causal colocalization and transcriptomic dysregulation converge.
#
# Shows the intersection of the SuSiE-COLOC gene set (368 genes,
# gene_level_coloc.csv coloc_best_susie_pp4 > 0.5) with the bulk Tier-1 DEGs
# (1,853 genes, canonical_deg_results.csv padj<0.05 & |logFC|>0.5).
#
#   x = bulk shrunk_logFC (ashr, canonical effect-size axis)
#   y = COLOC PP.H4 (SuSiE-best; axis labeled "COLOC PP.H4" because the
#       case-study set also includes abf-only genes)
#
# Convergent genes (SuSiE-COLOC AND Tier-1 DEG) are colored by direction
# (up = magenta, down = blue); other SuSiE-COLOC genes are grey. Case-study
# genes are directly labeled.
#
# CAVEAT (baked into legend): SuSiE is the canonical posterior (368 genes at
# PP.H4>0.5); the coloc.abf fallback calls 618. Axis is "COLOC PP.H4", NOT
# "SuSiE-", per the ledger.
#
# ALL numbers come from canonical files named in the Fig 4 number ledger:
#   GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv
#   RNA-seq/.../integration/results/integration/canonical_deg_results.csv
# NEVER from atlas *_coloc_pp4 convenience columns.
#
# Output: figures/main/fig4_validation/coloc_deg_convergence.pdf
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ===========================================================================
# 1. Load canonical files (ledger-named; no convenience columns)
# ===========================================================================
coloc_f <- file.path(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
deg_f   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")

coloc <- fread(coloc_f)
deg   <- fread(deg_f)

# ===========================================================================
# 2. Define the two evidence sets (ledger-verified counts)
# ===========================================================================
SUSIE_CUT <- 0.5   # SuSiE canonical PP.H4 threshold (ledger: 368 genes)
LFC_CUT   <- 0.5   # Tier-1 |logFC| floor

# Caveat counts over the FULL table (matches the ledger exactly; one abf-only
# hit carries an empty symbol). SuSiE = 368 canonical; abf fallback = 618.
n_susie <- sum(coloc$coloc_best_susie_pp4 > SUSIE_CUT, na.rm = TRUE)   # 368
n_abf   <- sum(coloc$coloc_best_pp4       > SUSIE_CUT, na.rm = TRUE)   # 618 (fallback)

# gene_level_coloc.csv `gene` column IS the HGNC symbol; drop the empty-symbol row
coloc <- coloc[gene != "" & !is.na(gene)]
# Tier-1 on the ASHR canonical axis (lfsr<0.05 & |shrunk_logFC|>0.5 = 1,433),
# to MATCH the plotted shrunk_logFC x-axis + the +/-0.5 guides. (Raw padj-based
# Tier-1 = 1,853, but mixing raw-LFC membership with a shrunk axis put 4
# "convergent" points inside the +/-0.5 band — verifier flag — so use ashr here.)
deg[, is_tier1 := lfsr < 0.05 & abs(shrunk_logFC) > LFC_CUT]
n_deg      <- sum(deg$is_tier1, na.rm = TRUE)                              # 1,433 (ashr)
n_deg_up   <- sum(deg$is_tier1 & deg$shrunk_logFC > 0, na.rm = TRUE)
n_deg_down <- sum(deg$is_tier1 & deg$shrunk_logFC < 0, na.rm = TRUE)

# A handful of HGNC symbols map to >1 Ensembl row in the canonical DEG table
# (e.g. CFHR5, HLA-DPB1) — collapse to one row per symbol BEFORE the COLOC join
# so each gene plots once. Prefer the most-significant row (lowest lfsr),
# tie-broken by the strongest |shrunk_logFC|. This keeps the convergent count
# at 18 and picks the Tier-1 row when any duplicate qualifies.
deg[, abs_slfc := abs(shrunk_logFC)]
setorder(deg, lfsr, -abs_slfc)
deg <- deg[!duplicated(symbol)]

# ===========================================================================
# 3. Build plot table: SuSiE-COLOC genes joined to canonical DEG by symbol
# ===========================================================================
susie <- coloc[coloc_best_susie_pp4 > SUSIE_CUT,
               .(symbol = gene,
                 coloc_pp4 = coloc_best_susie_pp4,
                 best_gwas = coloc_best_susie_gwas)]

dat <- merge(susie,
             deg[, .(symbol, shrunk_logFC, logFC, padj, lfsr, is_tier1)],
             by = "symbol", all.x = TRUE)

# Convergent = SuSiE-COLOC AND Tier-1 DEG (the headline intersection)
dat[, convergent := !is.na(is_tier1) & is_tier1]
dat[, status := fcase(
  convergent & shrunk_logFC > 0, "Convergent (up)",
  convergent & shrunk_logFC < 0, "Convergent (down)",
  default = "COLOC only"
)]

# Effect-size axis: prefer the canonical ashr shrunk_logFC; genes never tested
# in the DEG table (no canonical effect size) are dropped from the scatter.
dat <- dat[!is.na(shrunk_logFC)]

n_convergent <- sum(dat$convergent)
n_conv_up    <- sum(dat$status == "Convergent (up)")
n_conv_down  <- sum(dat$status == "Convergent (down)")

# ===========================================================================
# 4. Label set: all convergent genes + sub-threshold case studies
# ===========================================================================
# Convergent genes worth naming (top by PP.H4 * |effect|, plus key case studies)
case_convergent <- c("HKDC1", "FADS2", "SORT1", "CETP", "CFHR5", "EFHD1")
dat[, rank_score := coloc_pp4 * abs(shrunk_logFC)]
setorder(dat, -rank_score)
top_conv <- head(dat[convergent == TRUE, symbol], 12)
label_genes <- unique(c(case_convergent, top_conv))
dat[, show_label := symbol %in% label_genes]

# RORA is a strong COLOC hit (SuSiE 0.998) but a sub-threshold DEG
# (shrunk_logFC -0.368, below the 0.5 floor) — annotate it as a distinct marker
# to illustrate "genetics present, transcription below the Tier-1 cutoff".
rora <- dat[symbol == "RORA"]
dat[symbol == "RORA", show_label := TRUE]

# ===========================================================================
# 5. Scatter
# ===========================================================================
status_cols <- c("Convergent (up)"   = unname(masld_colors$up),    # #C9265E magenta
                 "Convergent (down)" = unname(masld_colors$down),  # #1565C0 blue
                 "COLOC only"        = unname(masld_colors$ns))    # #9E9E9E grey

# Plot grey background first, convergent on top
bg <- dat[convergent == FALSE]
fg <- dat[convergent == TRUE]

x_rng <- range(dat$shrunk_logFC)
x_pad <- diff(x_rng) * 0.10

p <- ggplot() +
  # threshold guides
  geom_hline(yintercept = SUSIE_CUT, linetype = "dashed",
             color = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed",
             color = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.25) +
  # COLOC-only (grey) points
  geom_point(data = bg, aes(x = shrunk_logFC, y = coloc_pp4),
             color = status_cols["COLOC only"], size = 0.7, alpha = 0.5) +
  # RORA: distinct open marker — strong genetics, sub-threshold transcription
  geom_point(data = rora, aes(x = shrunk_logFC, y = coloc_pp4),
             shape = 21, fill = NA, color = "grey25", stroke = 0.5, size = 2.1) +
  # convergent points (direction-colored), outlined
  geom_point(data = fg, aes(x = shrunk_logFC, y = coloc_pp4, fill = status),
             shape = 21, color = "white", stroke = 0.25, size = 2.1, alpha = 0.95) +
  scale_fill_manual(values = status_cols, name = NULL,
                    breaks = c("Convergent (up)", "Convergent (down)")) +
  geom_text_repel(
    data = dat[show_label == TRUE & symbol != "RORA"],
    aes(x = shrunk_logFC, y = coloc_pp4, label = symbol),
    size = 1.9, fontface = "italic", color = "grey15",
    max.overlaps = 40, segment.size = 0.2, segment.color = "grey60",
    min.segment.length = 0, box.padding = 0.4, point.padding = 0.25,
    force = 5, force_pull = 0.4, seed = 7
  ) +
  # RORA labeled separately so its sub-threshold caveat reads clearly
  geom_text_repel(
    data = rora, aes(x = shrunk_logFC, y = coloc_pp4,
                     label = "RORA (sub-threshold)"),
    size = 1.8, fontface = "italic", color = "grey30",
    nudge_x = 0.12, nudge_y = -0.085, segment.size = 0.2,
    segment.color = "grey55", min.segment.length = 0,
    box.padding = 0.3, seed = 7
  ) +
  scale_y_continuous(limits = c(0.4, 1.02),
                     breaks = c(0.5, 0.75, 1.0),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_x_continuous(limits = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Bulk log2FC (disease vs control, ashr-shrunk)",
       y = "COLOC PP.H4") +
  theme_masld(base_size = 7) + theme_pub() +
  theme(legend.position = c(0.015, 0.02),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.spacing.y = unit(0.02, "cm"),
        plot.margin = margin(4, 6, 3, 4))

# Convergence count annotation (Euler-style summary) — placed in the open
# lower-center band, clear of the labeled convergent points at the top.
conv_lab <- sprintf("SuSiE-COLOC (%d)  ∩  Tier-1 DEG (%d, ashr)  =  %d convergent",
                    n_susie, n_deg, n_convergent)
p <- p +
  annotate("text", x = 0, y = 0.59, label = conv_lab,
           hjust = 0.5, vjust = 1, size = 2.1, color = "grey20")

# ===========================================================================
# 6. Save
# ===========================================================================
out_pdf <- file.path(FIG4_DIR, "_supp", "coloc_deg_convergence.pdf")
dir.create(FIG4_DIR, recursive = TRUE, showWarnings = FALSE)
ggsave(out_pdf, p, width = 90, height = 78, units = "mm", device = cairo_pdf)

# ===========================================================================
# 7. Stats / legend text to stdout (NOT printed on the panel)
# ===========================================================================
message("================ coloc_deg_convergence.R ================")
message(sprintf("SuSiE-COLOC genes (gene_level_coloc.csv, coloc_best_susie_pp4>%.1f): %d",
                SUSIE_CUT, n_susie))
message(sprintf("  CAVEAT: SuSiE is canonical (%d genes); coloc.abf fallback calls %d at PP.H4>%.1f.",
                n_susie, n_abf, SUSIE_CUT))
message(sprintf("Bulk Tier-1 DEGs (canonical_deg_results.csv, padj<0.05 & |logFC|>%.1f): %d (%d up / %d down)",
                LFC_CUT, n_deg, n_deg_up, n_deg_down))
message(sprintf("CONVERGENT (SuSiE-COLOC AND Tier-1 DEG): %d genes (%d up / %d down)",
                n_convergent, n_conv_up, n_conv_down))
message("Axis label is 'COLOC PP.H4' (NOT 'SuSiE-') because the case-study set")
message("  includes abf-only genes (e.g. THRB abf 1.000 with NO SuSiE).")
message("Effect-size axis = ashr-shrunk_logFC (canonical) from canonical_deg_results.csv.")
message("\nConvergent genes (PP.H4 [SuSiE], shrunk_logFC, best GWAS):")
setorder(dat, -coloc_pp4)
for (i in which(dat$convergent)) {
  message(sprintf("  %-9s  PP.H4=%.3f  log2FC=%+.3f  (%s)",
                  dat$symbol[i], dat$coloc_pp4[i], dat$shrunk_logFC[i], dat$best_gwas[i]))
}
message(sprintf("\nNotable sub-threshold COLOC case study: RORA SuSiE PP.H4=%.3f but shrunk_logFC=%.3f (below |%.1f| Tier-1 floor; not convergent).",
                coloc[gene == "RORA", coloc_best_susie_pp4],
                deg[symbol == "RORA", shrunk_logFC], LFC_CUT))
message(sprintf("Saved: %s", out_pdf))
message("=========================================================")
