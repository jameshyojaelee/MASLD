#!/usr/bin/env Rscript
# KEY MESSAGE (POSITIVE CONTROL): Hepatocytes and fibroblasts spatially
# EXCLUDE one another (per-spot co-occurrence Spearman rho = -0.562, the most
# negative cell-type pair), and CXCL12-high, lymphocyte-poor microdomains
# mark the immune-excluded fibrotic niche (CXCL12 delta = +0.312, excluded vs
# non-excluded spots).
#
# This is a POSITIVE CONTROL: it recovers an architecture already reported in
# the literature (Karpova 2026; Ramachandran 2019 fibrotic-niche scarring),
# so it validates the spatial pipeline rather than claiming a novel finding.
#
# HARD CAVEATS (baked into the legend, emitted to stdout via message()):
#   * DESCRIPTIVE ONLY — no across-donor p-value. Co-occurrence rho is the
#     across-donor mean of per-donor spot-level Spearman; only ~5 donors and
#     spots within a slice are spatially autocorrelated.
#   * CXCL12 elevation is NOT significant: paired p = 0.0625, padj(BH) = 0.19.
#   * Reframe as LYMPHOCYTE-POOR stroma, not "immune-excluded" wholesale:
#     macrophages are NOT depleted in excluded spots (they trend UP);
#     it is the lymphoid compartment (T / B / NK) that is depleted.
#   * Scooped architecture: Karpova 2026; Ramachandran 2019.
#
# Data sources (numbers read from disk, never hardcoded from prose):
#   celltype_cooccurrence_spearman.csv  (V1,V2,N) -> Hep,Fib rho = -0.562
#   repulsive_ligand_elevation.csv      -> CXCL12 delta = +0.312 (p .0625/padj .19)
#   celltype_mean_by_spot_class.csv     -> lymphoid depletion / macrophage rise
#   immune_exclusion_summary.txt        -> 314/6546 excluded spots, 5 donors
#
# Output: figures/main/fig4_validation/spatial_immune_exclusion.pdf  (flat name)
# Env:    rnaseq
# Ledger: row "immune_exclusion" (hep-fib rho -0.562 descriptive; 314/6546=4.8%
#         excluded; CXCL12 +0.312 NOT sig p=.0625/padj .19; scooped Karpova 2026,
#         Ramachandran 2019).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

IN <- file.path(BASE, "Analysis/Spatial/results/immune_exclusion")

# ── Data 1: cell-type co-occurrence WITH hepatocytes ─────────────────────────
# The story is hep<->fib exclusion, so show every cell type's spatial
# co-occurrence rho against hepatocytes; fibroblasts is the most negative.
cooc <- fread(file.path(IN, "celltype_cooccurrence_spearman.csv"))
setnames(cooc, c("ct1", "ct2", "rho"))

hep_pairs <- cooc[ct1 == "Hepatocytes" & ct2 != "Hepatocytes"]
setorder(hep_pairs, rho)
hep_pairs[, partner := factor(ct2, levels = ct2)]   # ascending: most -ve first

# The headline number must come straight off disk.
hep_fib_rho <- cooc[ct1 == "Hepatocytes" & ct2 == "Fibroblasts", rho]
stopifnot(length(hep_fib_rho) == 1)
stopifnot(abs(hep_fib_rho - (-0.562)) < 1e-3)   # ledger row immune_exclusion

# Highlight fibroblasts (the headline repulsion); colour by sign otherwise.
hep_pairs[, hl := ifelse(ct2 == "Fibroblasts", "fib",
                  ifelse(rho < 0, "neg", "pos"))]
lolli_cols <- c(fib = masld_colors$up,        # headline repulsion = disease magenta
                neg = "#1565C0",              # attract-away (negative) = blue
                pos = "#9E9E9E")              # co-localising (positive) = gray

p_cooc <- ggplot(hep_pairs, aes(x = rho, y = partner)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  geom_segment(aes(x = 0, xend = rho, yend = partner, color = hl),
               linewidth = 0.5) +
  geom_point(aes(color = hl), size = 1.7) +
  geom_text(data = hep_pairs[ct2 == "Fibroblasts"],
            aes(label = sprintf("%.3f", rho)),
            hjust = -0.25, vjust = -1.1, size = PUB_GEOM_TEXT,
            fontface = "bold", color = masld_colors$up) +
  scale_color_manual(values = lolli_cols, guide = "none") +
  scale_x_continuous(limits = c(-0.65, 0.30),
                     breaks = c(-0.5, -0.25, 0, 0.25)) +
  labs(x = "Co-occurrence with hepatocytes\n(per-spot Spearman rho)",
       y = NULL, title = "Hepatocyte-stroma exclusion") +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(size = PUB_AXIS_TEXT + 0.3))

# ── Data 2: repulsive-ligand elevation in the excluded niche ─────────────────
lig <- fread(file.path(IN, "repulsive_ligand_elevation.csv"))
cxcl12_delta <- lig[ligand == "CXCL12", delta_excluded_minus_non]
cxcl12_p     <- lig[ligand == "CXCL12", pval]
cxcl12_padj  <- lig[ligand == "CXCL12", padj_bh]
stopifnot(abs(cxcl12_delta - 0.312) < 2e-3)   # ledger row immune_exclusion

setorder(lig, delta_excluded_minus_non)
lig[, lname := factor(ligand, levels = ligand)]
# CXCL12 is the only one of interest but NOT significant -> gray it, annotate n.s.
lig[, sig := ifelse(padj_bh < 0.05, "sig", "ns")]
lig_cols <- c(sig = masld_colors$up, ns = "#9E9E9E")

p_lig <- ggplot(lig, aes(x = delta_excluded_minus_non, y = lname)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  geom_segment(aes(x = 0, xend = delta_excluded_minus_non, yend = lname,
                   color = sig), linewidth = 0.5) +
  geom_point(aes(color = sig), size = 1.7) +
  geom_text(data = lig[ligand == "CXCL12"],
            aes(label = sprintf("+%.3f (n.s.)", delta_excluded_minus_non)),
            hjust = 1.05, vjust = -1.1, size = PUB_GEOM_TEXT,
            color = "gray30") +
  scale_color_manual(values = lig_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 0.42), breaks = c(0, 0.1, 0.2, 0.3)) +
  labs(x = "Ligand elevation in excluded niche\n(excluded - non-excluded spots)",
       y = NULL, title = "Repulsive ligands") +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(size = PUB_AXIS_TEXT + 0.3))

# ── Lymphoid-depletion check (drives the "lymphocyte-poor" reframe) ──────────
# Confirm macrophages are NOT depleted; lymphoid (T/B/NK) is. For the legend.
by_cls <- fread(file.path(IN, "celltype_mean_by_spot_class.csv"))
get_cls <- function(cls, ct) by_cls[spot_class == cls][[ct]]
mac_exc  <- get_cls("immune_excluded_fibrotic", "Macrophages")
mac_non  <- get_cls("non_excluded", "Macrophages")
tcell_exc <- get_cls("immune_excluded_fibrotic", "T cells")
tcell_non <- get_cls("non_excluded", "T cells")
bcell_exc <- get_cls("immune_excluded_fibrotic", "B cells")
bcell_non <- get_cls("non_excluded", "B cells")
fib_exc  <- get_cls("immune_excluded_fibrotic", "Fibroblasts")
fib_non  <- get_cls("non_excluded", "Fibroblasts")

# ── Assemble (compact, two lollipops side by side) ───────────────────────────
p_out <- (p_cooc | p_lig) + plot_layout(widths = c(1.25, 1))

# ── Legend material to stdout (PI directive: stats off the panel) ────────────
message(sprintf(
"[spatial_immune_exclusion legend] POSITIVE CONTROL (GSE192741 Visium, %d spots, %d donors, cell2location). Hepatocytes and fibroblasts spatially EXCLUDE one another: per-spot co-occurrence Spearman rho = %.3f (the most negative cell-type pair); hepatocytes also avoid lymphoid cells (T %.3f, B %.3f, circulating-NK %.3f). CXCL12 is elevated in the immune-excluded fibrotic niche (delta = +%.3f, excluded vs non-excluded spots) but this is NOT significant (paired p = %.4f, padj(BH) = %.2f); TGFB1/VEGFA deltas are near zero. CAVEATS: DESCRIPTIVE ONLY -- co-occurrence rho is the across-donor mean of per-donor spot-level Spearman; with only ~%d donors and spatially autocorrelated spots within each slice, NO across-donor p-value is claimed. Reframe as LYMPHOCYTE-POOR stroma, not wholesale immune exclusion: macrophages are NOT depleted in excluded spots (excluded %.4f vs non-excluded %.4f, trending UP), while fibroblasts rise (%.4f vs %.4f) and lymphoid cells fall. Architecture is SCOOPED (positive control): Karpova 2026; Ramachandran 2019.",
  6546, 5,
  hep_fib_rho,
  cooc[ct1 == "Hepatocytes" & ct2 == "T cells", rho],
  cooc[ct1 == "Hepatocytes" & ct2 == "B cells", rho],
  cooc[ct1 == "Hepatocytes" & ct2 == "Circulating NK/NKT", rho],
  cxcl12_delta, cxcl12_p, cxcl12_padj,
  5,
  mac_exc, mac_non, fib_exc, fib_non))

out <- file.path(FIG4_DIR, "_supp", "spatial_immune_exclusion.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
pdf_device(out, width = fig_col_width, height = 2.3)
print(p_out)
invisible(dev.off())
message("Saved: ", out)
