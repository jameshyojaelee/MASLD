##############################################################################
# fig5_panel_calibration.R — Fig 5 PANEL 5d (REBUILT 2026-06-22)
#
# CALIBRATION of the multi-evidence CONVERGENCE SCORE (Script 46d) against
# INDEPENDENT external MASLD truth panels — vs CHANCE, not vs our old method.
#
# Objective: show the convergence ranking recovers diverse, independent external
# truth ABOVE CHANCE (so the ranking is real / non-circular). NOT a horse race.
# The previous "vs prior method (46b, retired)" comparison was dropped — beating
# our own deprecated method is a development artifact, not a result.
#
# Three external truth panels, each a different evidence type:
#   Govaere 2020   — published bulk MASLD expression signature   (expression)
#   NIDDK 2024     — curated MASLD clinical drug targets          (clinical)
#   Open Targets   — MASLD genetic associations (NASH EFO_1001249) (genetic)
#
# HONEST framing (caption, emitted to stdout — NOT on the panel):
#  - The convergence score recovers EXPRESSION (Govaere) and CLINICAL (NIDDK)
#    truth strongly, and PURELY-GENETIC truth (Open Targets) weakly/near-chance.
#    This is BY DESIGN: the score rewards genes with CONVERGENT multi-modal
#    support, so a gene with only a genetic association (no expression/protein/
#    single-cell signal) scores low. The genetic panel tests what the score is
#    deliberately NOT built to do — it characterizes the score, it is not a flaw.
#  - Open Targets panel uses NASH EFO_1001249 only; the EFO_0004612 set is HDL
#    (not NAFLD; see memory) and is dropped.
#  - Govaere shares samples with the bulk layer (GSE135251); a held-out
#    sensitivity (bulk DEG recomputed excluding Govaere) moves the bulk-layer
#    AUROC only 0.905 -> 0.891, so the recovery is not driven by sample overlap.
#  - Govaere / NIDDK / Open Targets are TRUTH sets being predicted, NOT rival
#    tools. NEVER imply the score "beats" those external resources.
##############################################################################

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANDIR <- file.path(FIG5_DIR, "panels"); dir.create(PANDIR, recursive = TRUE, showWarnings = FALSE)

# ── Convergence score + external truth panels ────────────────────────────────
A  <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
ce <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/convergence_evidence.csv"))
A  <- merge(A, ce[, .(human_symbol, conv = convergence_score)], by = "human_symbol", all.x = TRUE)
A  <- A[is.finite(conv)]                                    # scored universe (matches 46d benchmark)
rd <- function(p) { L <- readLines(p); L <- L[!grepl("^#", L)]
                    unique(fread(text = paste(L, collapse = "\n"))$gene_symbol) }
ot <- { L <- readLines(file.path(BASE, "data/published_gene_panels/opentargets_masld_2025.tsv"))
        d <- fread(text = paste(L[!grepl("^#", L)], collapse = "\n"))
        unique(d[source_disease == "NASH_EFO_1001249", gene_symbol]) }   # drop HDL EFO_0004612
panels <- list(
  Govaere = list(genes = rd(file.path(BASE, "data/published_gene_panels/govaere_2020_panel.tsv")),
                 label = "Govaere signature\n(expression)"),
  NIDDK   = list(genes = rd(file.path(BASE, "data/published_gene_panels/niddk_pipeline_2024.tsv")),
                 label = "NIDDK targets\n(clinical)"),
  OpenTargets = list(genes = ot, label = "Open Targets\n(genetic)"))

aurocf <- function(s, lab) { r <- rank(s, ties.method = "average")
  n1 <- sum(lab); n0 <- sum(!lab); (sum(r[lab]) - n1*(n1+1)/2) / (n1*n0) }
# Hanley & McNeil (1982) analytic 95% CI — the standard AUROC SE; correctly
# wide for small positive sets (a stratified bootstrap under-states it because
# the few panel genes are consistently top-ranked).
hanley <- function(a, n1, n0) { Q1 <- a/(2-a); Q2 <- 2*a^2/(1+a)
  se <- sqrt((a*(1-a) + (n1-1)*(Q1-a^2) + (n0-1)*(Q2-a^2)) / (n1*n0))
  c(lo = max(0, a - 1.96*se), hi = min(1, a + 1.96*se)) }
res <- rbindlist(lapply(names(panels), function(pn) {
  lab <- A$human_symbol %in% panels[[pn]]$genes; s <- A$conv
  au  <- aurocf(s, lab); ci <- hanley(au, sum(lab), sum(!lab))
  p   <- suppressWarnings(wilcox.test(s[lab], s[!lab], alternative = "greater")$p.value)
  data.table(panel = pn, label = panels[[pn]]$label, n_pos = sum(lab),
             auroc = au, ci_lo = ci["lo"], ci_hi = ci["hi"], p = p) }))
res[, star := fifelse(p < 1e-3, "***", fifelse(p < 1e-2, "**", fifelse(p < 0.05, "*", "ns")))]
res[, label := factor(label, levels = res[order(-auroc), label])]
fwrite(res, file.path(PANDIR, "fig5d_calibration_source.csv"))

# ── Panel 5d — convergence-score AUROC ± bootstrap CI per panel, vs chance ────
p5d <- ggplot(res, aes(x = label, y = auroc)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "grey55") +
  geom_col(width = 0.62, fill = masld_colors$mash, alpha = 0.95) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.18, linewidth = 0.4, color = "black") +
  geom_text(aes(label = sprintf("%.2f%s", auroc, star), y = ci_hi), vjust = -0.5, size = 2.3, color = "black") +
  annotate("text", x = 0.55, y = 0.515, label = "chance", hjust = 0, vjust = 0, size = 2.0, color = "grey45") +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "AUROC (external-panel recovery)", title = "d") +
  theme_masld(base_size = 7) +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.text.x = element_text(size = 6.2, lineheight = 0.85))
ggsave(file.path(PANDIR, "fig5d_calibration.pdf"), p5d, width = 3.4, height = 3.2, device = cairo_pdf)

cat("Saved panels/fig5d_calibration.pdf + source.csv\n\n")
print(res[, .(panel, n_pos, auroc = round(auroc,3), ci = sprintf("[%.2f,%.2f]", ci_lo, ci_hi), p = signif(p,2), star)])
cat("\nCAPTION (not on panel): The convergence score recovers independent external MASLD truth above\n",
    "chance for expression (Govaere) and clinical-target (NIDDK) panels; recovery of the purely-genetic\n",
    "Open Targets panel is near chance BY DESIGN (the score rewards multi-modal convergence, so genes with\n",
    "only a genetic association score low). Panels are held-out TRUTH, not competing tools.\n")
