#!/usr/bin/env Rscript
# ============================================================================
# figS06_vacca_benchmark.R — our cross-species results vs Vacca 2024 (LITMUS)
#
# Panel A (model-level VALIDATION): LITMUS's DHPS is two-armed (metabolic +
#   fibrotic). A pathway-level PEP reproduces the FIBROTIC arm across LITMUS's
#   own 33 mouse models (Spearman 0.62 full-pathway / 0.54 on the 29 pathways we
#   share; metabolic arm does NOT reproduce, 0.29 NS — deferred to LITMUS). On
#   that validated fibrotic axis our 4 diets land sensibly: FPC (Western) is the
#   most human-proximal, reproducing LITMUS's "Western closest" result; MCD is
#   genuinely fibrosis-strong (its metabolic demerit is the un-reproduced arm).
# Panel B (gene-level CONTRAST-MATCHED concordance): our human signatures ×
#   LITMUS's stage-partitioned 951 signature. The progression contrast (NAFL-vs-
#   NASH) matches LITMUS's progression bin best (OR 31.7), disease-vs-control
#   matches the early bin only weakly (OR 1.7, near-null specificity). A
#   threshold-free, non-circular anchor (genes concordant in ≥1 mouse diet) is
#   4.3× enriched in LITMUS-951 (p=1e-72), robust to Govaere holdout + shuffle.
#
# ⚠ Gene-level concordance shares human input (LITMUS EPoS = our Govaere cohort);
#   we drop all human-vs-human circular metrics (RRHO / signed Spearman) and report
#   the contrast-matched OR + the cross-species (mouse-anchored) non-circular leg.
# Data: Analysis/Cross_Species_Concordance/results/vacca_benchmark/ (scripts 11, 12).
# Cite: Vacca et al. 2024 Nat Metab, DOI 10.1038/s42255-024-01043-6 (CC BY).
# Output: figures/supplementary/figS06_cross_species/figS06_vacca_benchmark.pdf
# Env: rnaseq
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

VB   <- file.path(BASE, "Analysis/Cross_Species_Concordance/results/vacca_benchmark")
MAG  <- "#C9265E"; TEAL <- "#00695C"; GREY <- "#9E9E9E"; BLUE <- "#1565C0"

# ── Panel A: fibrotic-arm validation + our diets on the validated axis ────────
# 33 LITMUS models: x = our reproduced fibrotic-PEP (shared-29 axis, the one our
# diets live on), y = LITMUS published fibrotic-arm DHPS. Our 4 diets carry no
# published DHPS -> shown as a labelled rug + bootstrap-CI on the SAME x-axis.
litm  <- fread(file.path(VB, "q1_litmus_models_fibrotic_pep_shared.csv"))   # model,fibro_pep,diet_group,f_DHPS
arms  <- fread(file.path(VB, "q1_fibrotic_validation_arms.csv"))
our   <- fread(file.path(VB, "q1_our_diets_fibrotic_pep.csv"))
rho_shared <- cor(litm$fibro_pep, litm$f_DHPS, method = "spearman", use = "complete.obs")
rho_full   <- arms[grepl("^MATCHED: fibrotic", test)]$spearman
p_full     <- arms[grepl("^MATCHED: fibrotic", test)]$p_value

west_grp <- function(g) grepl("WD|GAN|AMLN|AMLD", g)
litm[, west := west_grp(diet_group)]
litm_med <- median(litm$fibro_pep, na.rm = TRUE)
our[, lab := our_diet]
our[, dy  := factor(our_diet, levels = our[order(fibro_pep)]$our_diet)]    # FPC at top
xlim <- c(min(c(litm$fibro_pep, our$fibro_ci_lo)) - 0.02, 1.0)
xbreaks <- c(0.3, 0.5, 0.7, 0.9)

# a_top: our 4 diets, each its own row, point + bootstrap 95% CI; shares x with a_bot
a_top <- ggplot(our) +
  geom_vline(xintercept = litm_med, linetype = "22", colour = "grey70", linewidth = 0.35) +
  geom_segment(aes(x = fibro_ci_lo, xend = fibro_ci_hi, y = dy, yend = dy),
               colour = MAG, linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(x = fibro_pep, y = dy), colour = MAG, size = 2.4) +
  scale_x_continuous(limits = xlim, breaks = xbreaks) +
  labs(tag = "a", x = NULL, y = NULL, title = "our 4 diets (no published DHPS)") +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_text(colour = "black", face = "bold.italic"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(size = PUB_GEOM_TEXT * 2.6, colour = "grey35", hjust = 0),
        plot.tag = element_text(size = 11, face = "bold"))

# a_bot: validation scatter — reproduced fibrotic-PEP vs LITMUS published f_DHPS
a_bot <- ggplot(litm) +
  geom_vline(xintercept = litm_med, linetype = "22", colour = "grey70", linewidth = 0.35) +
  geom_smooth(aes(x = fibro_pep, y = f_DHPS), method = "lm", se = FALSE,
              colour = "grey75", linewidth = 0.4, linetype = "22") +
  geom_point(aes(x = fibro_pep, y = f_DHPS, shape = west),
             colour = GREY, fill = "grey80", size = 1.6, alpha = 0.9, stroke = 0.3) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
  annotate("text", x = 1.0, y = 0.04, hjust = 1, vjust = 0, size = PUB_GEOM_TEXT,
           colour = "black", label = sprintf("LITMUS models: ρ = %.2f", rho_shared)) +
  scale_x_continuous(limits = xlim, breaks = xbreaks) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
  labs(x = "fibrotic pathway-proximity (PEP corr. to human Severe-vs-Mild)",
       y = "LITMUS published\nfibrotic-arm DHPS") +
  theme_masld() + theme_pub() +
  theme(axis.text = element_text(colour = "black"))

p_a <- a_top / a_bot + plot_layout(heights = c(0.62, 1.5))

# ── Panel B: contrast-matched enrichment matrix (Fisher OR) ───────────────────
mat <- fread(file.path(VB, "q2_contrast_matrix.csv"))
mat <- mat[human_sig %in% c("H2_naflnash", "H1_dvc")]                     # headline rows (current method)
rowlab <- c(H2_naflnash = "NAFL→NASH\n(progression)", H1_dvc = "disease\nvs control")
collab <- c(V_early_271 = "early\n(271)", V_prog_526 = "progression\n(526)", V_all_951 = "all\n(951)")
mat[, rl := factor(rowlab[human_sig], levels = c(rowlab["H1_dvc"], rowlab["H2_naflnash"]))]
mat[, cl := factor(collab[as.character(vacca_target)],
                   levels = c(collab["V_early_271"], collab["V_prog_526"], collab["V_all_951"]))]
mat[, l2or := log2(OR)]
mat[, txt  := ifelse(OR >= 10, sprintf("%.0f×", OR), sprintf("%.1f×", OR))]

p_b <- ggplot(mat, aes(x = cl, y = rl, fill = l2or)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = txt), colour = "black", fontface = "bold", size = PUB_GEOM_TEXT) +
  scale_fill_gradient2(low = "#B33018", mid = "white", high = TEAL, midpoint = 0,
                       name = expression(log[2]*"(OR)"),
                       breaks = c(0, 2, 4), limits = c(-1.2, max(mat$l2or) + 0.2)) +
  scale_x_discrete(position = "top") +
  labs(tag = "b", x = "LITMUS-951 stage bin", y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text = element_text(colour = "black"),
        axis.text.x.top = element_text(colour = "black"),
        legend.position = "right", legend.key.width = unit(0.28, "cm"),
        plot.tag = element_text(size = 11, face = "bold"),
        panel.grid = element_blank())

# ── Caption / provenance → stdout (stats live here, not on the panels) ────────
anc  <- fread(file.path(VB, "q2_noncircular_anchor.csv"))[metric == "fisher_nconc1_vs_litmus951"]
gov  <- fread(file.path(VB, "q2_govaere_holdout.csv"))
nul  <- fread(file.path(VB, "q2_shuffle_null.csv"))
gprog<- dcast(gov[vacca_target == "V_prog_526"], vacca_target ~ sig_set, value.var = "OR")
diag <- mat[human_sig == "H2_naflnash" & vacca_target == "V_prog_526"]
ctrl <- mat[human_sig == "H1_dvc" & vacca_target == "V_early_271"]
fpc  <- our[lab == "FPC"]; mcd <- our[lab == "MCD"]

message(strrep("=", 78))
message("FIGURE LEGEND (figS06_vacca_benchmark)")
message(strrep("=", 78))
message(sprintf(
"Benchmarking our cross-species analysis against the LITMUS murine-model ranking
(Vacca et al. 2024 Nat Metab, DOI 10.1038/s42255-024-01043-6; CC BY). LITMUS's DHPS
is a two-armed KEGG-PATHWAY construct. (a) A pathway-level PEP reproduces LITMUS's
FIBROTIC-arm DHPS across their own %d mouse models (Spearman %.2f over the %d shared
KEGG pathways; %.2f over the full pathway set, p=%.0e); the METABOLIC arm does not
reproduce (0.29, NS) and is deferred to LITMUS. On this validated fibrotic axis our 4
diets land sensibly: FPC (Western) is the most human-proximal (PEP %.2f, bootstrap 95%%
CI [%.2f, %.2f]; 0 of %d LITMUS models above), reproducing LITMUS's 'Western diets
closest to human'. MCD is genuinely fibrosis-strong (fibrotic rank %d; LITMUS MCD/CDD
f_DHPS~0.75) — its low metabolic standing is the un-reproduced arm, not an error. (The
FPC point estimate is leverage-sensitive — few-pathway-driven — hence the wide CI; we
present consistency with, and defer the precise ranking to, LITMUS.)
(b) Contrast-matched concordance of our human signatures with LITMUS's stage-partitioned
951 signature: the progression contrast (NAFL→NASH) matches LITMUS's progression bin best
(OR %.1f, p=%.0e), while disease-vs-control matches the early bin only weakly (OR %.1f) —
apples-to-apples lifts the overlap, so the modest raw 172-gene overlap was a contrast
mismatch, not a ceiling. A threshold-free, NON-CIRCULAR anchor (genes concordant in ≥1
mouse diet, no human cutoff) is %.1f× enriched in LITMUS-951 (Fisher p=%.0e, n=%d),
robust to dropping the shared Govaere cohort (progression OR %.1f→%.1f) and to 1,000×
label shuffle (p<0.001). Circular human-vs-human metrics (RRHO, signed Spearman) are
NOT reported.",
  nrow(litm), rho_shared, our[1]$n_shared, rho_full, p_full,
  fpc$fibro_pep, fpc$fibro_ci_lo, fpc$fibro_ci_hi, nrow(litm), mcd$fibro_rank,
  diag$OR, diag$fisher_p, ctrl$OR,
  anc$OR, anc$p, anc$n_both, gprog$FULL, gprog$NOGOV))
message("CAVEAT: gene-level (b) shares human input (LITMUS EPoS = our Govaere); the non-circular leg is mouse-anchored.")
message(strrep("=", 78))

# ── Assemble + save ──────────────────────────────────────────────────────────
p_out <- (p_a | p_b) + plot_layout(widths = c(1.55, 1.5))
out <- file.path(FIGS06_DIR, "figS06_vacca_benchmark.pdf")
dir.create(FIGS06_DIR, showWarnings = FALSE, recursive = TRUE)
save_fig(p_out, out, width = fig_full_width * 0.95, height = 2.5)
message("Saved: ", out)

# ── Supp tables → figS06 data/ ────────────────────────────────────────────────
dd <- file.path(FIGS06_DIR, "data"); dir.create(dd, showWarnings = FALSE, recursive = TRUE)
for (f in c("q1_fibrotic_validation_arms.csv", "q1_our_diets_fibrotic_pep.csv",
            "q1_fibrotic_per_model.csv", "q2_contrast_matrix.csv",
            "q2_noncircular_anchor.csv", "q2_govaere_holdout.csv"))
  file.copy(file.path(VB, f), file.path(dd, paste0("figS06_", f)), overwrite = TRUE)
message("Supp tables copied to: ", dd)
