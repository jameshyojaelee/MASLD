#!/usr/bin/env Rscript
# =============================================================================
# 83_phase5_figures.R — Phase-5 seqfunc SUPPLEMENTARY figures (MASLD house style)
# =============================================================================
# ADDITIVE ONLY. Builds supplementary panels from the SETTLED Phase-5 results
# (direction hardening + coding GoF/LoF benchmark). Does NOT touch 60-82, 46d,
# 78, 27a, or any convergence channel — the direction result is a bounded
# SUPPLEMENTARY confidence-recalibration finding, never a scored channel.
#
# Settled sources (locked; SF = GWAS/finemapping/results/seqfunc):
#   direction : SF/direction_hardened_results.csv + direction_hardened_verdict.json
#   coding    : SF/coding_hardening_v2.tsv + gof_lof_benchmark.json
#
# Panels rendered:
#   (1) direction_benchmark.pdf  — supervised model vs 4 zero-shot/positional
#        baselines, dashed 0.5 chance line. KEY MESSAGE: supervised meta-learning
#        does NOT beat zero-shot or a trivial positional baseline (thesis-reinforcing:
#        sequence-to-function models stay magnitude/mechanism-only). Method NOT novel (ExPecto,
#        Zhou 2018 Nat Genet).
#   (2) leakage_gap.pdf          — random-split 0.567 vs chromosome-held-out (LOCO)
#        0.507; the 0.059 gap = LD leakage, not biology.
#   (3) coding_gof_lof.pdf       — 42 coding effectors, per-model calls
#        (AlphaMissense / ESM1b / popEVE severity + LoGoFunc GoF/LoF/Neutral +
#        curated), spotlighting 3 canonical lipid effectors (PNPLA3 I148M /
#        TM6SF2 E167K / APOE C130R). Only PNPLA3 I148M is a true GoF neomorph
#        (TM6SF2 E167K = curated LoF, APOE C130R = isoform-functional); no
#        predictor emits a GoF call. Held-out GoF-vs-LoF auROC 0.905
#        [0.879, 0.930] annotated.
#
# HOOKS (clearly commented, below) left for:
#   - caQTL direction panel  (ChromBPNet / Borzoi caQTL arm — features landing)
#   - Decima cell-type panel (decima-celltype agent output — landing)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

SF      <- file.path(BASE, "GWAS/finemapping/results/seqfunc")
OUT_DIR <- file.path(BASE, "figures/supplementary/seqfunc")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

rendered <- character(0)

# Phase-6 (2026-07-11): the direction + Decima panels below are SUPERSEDED by
# src/88_phase6_figures.R (ToS-clean zero-shot). Guard save_fig so re-running this
# script cannot regenerate the retired PDFs or overwrite the new Decima panel.
# Only coding_gof_lof.pdf (Panel 3) is still owned here.
.RETIRED_PDFS <- c("direction_benchmark.pdf", "leakage_gap.pdf",
                   "qtl_modality_direction.pdf", "decima_celltype.pdf")
.orig_save_fig <- save_fig
save_fig <- function(p, path, ...) {
  if (basename(path) %in% .RETIRED_PDFS) {
    message("[SKIP retired -> see src/88] ", basename(path)); return(invisible(NULL))
  }
  .orig_save_fig(p, path, ...)
}

# =============================================================================
# PANEL 1 — direction_benchmark.pdf
# =============================================================================
# Forest / point-range: the tested supervised elastic-net (CERTIFIED_PRIMARY,
# distinct magenta) against the 4 pre-registered zero-shot / positional baselines
# (control gray #9E9E9E), all vs a dashed 0.5 chance line. CI whiskers from the
# block-bootstrap. Supervised sits AT chance, BELOW the zero-shot baselines and
# CI-overlapping the trivial positional baseline -> beats NONE (CI-separated).
# -----------------------------------------------------------------------------
dir_csv <- fread(file.path(SF, "direction_hardened_results.csv"))
# The caQTL re-run added caqtl/sqtl rows (multi-modality). Panels 1 & 2 are the
# CERTIFIED_PRIMARY eQTL expression benchmark only -> scope to that panel.
if ("panel" %in% names(dir_csv)) dir_csv <- dir_csv[panel == "eqtl_broadaway"]

# The tested supervised model + the 4 baselines named in the pre-registration.
# (agg_signed_mean and gbt_fixed are excluded from the panel: agg_signed_mean is
#  a redundant 5th baseline and gbt is a sensitivity arm; both are in the CSV.)
keep <- c(
  "elastic_net_nestedLOCO",          # supervised (the tested thing)
  "aggregated_signed_logSED",        # zero-shot baseline
  "fold_unanimous_rule",             # zero-shot baseline
  "logSED_gated_sign",               # zero-shot baseline
  "MAF_TSSdist_positional_logistic"  # trivial positional baseline
)
label_map <- c(
  elastic_net_nestedLOCO          = "Supervised elastic-net (LOCO)",
  aggregated_signed_logSED        = "Aggregated signed logSED",
  fold_unanimous_rule             = "Fold-unanimous sign rule",
  logSED_gated_sign               = "|logSED|-gated sign",
  MAF_TSSdist_positional_logistic = "MAF + TSS-distance (positional)"
)

d1 <- dir_csv[model %in% keep]
d1[, label   := label_map[model]]
d1[, is_supervised := model == "elastic_net_nestedLOCO"]
d1[, grp     := ifelse(is_supervised, "Supervised (tested)", "Zero-shot / positional baseline")]
# Order: best auROC at top, chance line reference; supervised falls mid-pack.
d1 <- d1[order(auroc)]
d1[, label := factor(label, levels = label)]

p1 <- ggplot(d1, aes(x = auroc, y = label, color = grp)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3,
             color = "#9E9E9E") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.35) +
  geom_point(aes(shape = grp), size = 1.7) +
  scale_color_manual(values = c("Supervised (tested)" = masld_colors$up,
                                "Zero-shot / positional baseline" = "#9E9E9E"),
                     name = NULL) +
  scale_shape_manual(values = c("Supervised (tested)" = 18,
                                "Zero-shot / positional baseline" = 16),
                     name = NULL) +
  scale_x_continuous(limits = c(0.40, 0.66), breaks = seq(0.40, 0.65, 0.05)) +
  labs(x = "Sign-concordance auROC (pooled OOF)", y = NULL,
       title = "Direction prediction: supervised vs zero-shot") +
  theme_masld() +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -2))

# CAPTION -> message() per house style (feedback-no-colored-fonts): no subtitle.
message("[direction_benchmark] CAPTION: Signed-direction prediction for MASLD ",
        "eQTL effects (hardened Broadaway panel, n=211, pre-registered NEGATIVE). ",
        "Points = block-bootstrap auROC, whiskers = 95% CI, dashed line = 0.5 chance. ",
        "The supervised elastic-net (nested leave-one-chromosome-out; magenta) sits ",
        "at chance (0.507 [0.437,0.575], perm-p=0.39) and CI-overlaps every zero-shot ",
        "and the trivial positional baseline -- it beats NONE (CI-separated). ",
        "Sequence-to-function model direction adds no signal beyond magnitude; method is not ",
        "novel (ExPecto, Zhou 2018 Nat Genet). Supplementary confidence-recalibration ",
        "only; APPLY-ONLY firewall -- never a scored convergence channel.")

save_fig(p1, file.path(OUT_DIR, "direction_benchmark.pdf"),
         width = 3.6, height = 2.4)
rendered <- c(rendered, "direction_benchmark.pdf")

# =============================================================================
# PANEL 2 — leakage_gap.pdf
# =============================================================================
# Compact 2-point comparison: random-split auROC (leaky) vs chromosome-held-out
# (LOCO, honest). The 0.059 drop is LD leakage, not biology.
# -----------------------------------------------------------------------------
prim <- dir_csv[model == "elastic_net_nestedLOCO"][1]
d2 <- data.table(
  split = factor(c("Random split\n(leaky)", "Chromosome-held-out\n(LOCO, honest)"),
                 levels = c("Random split\n(leaky)",
                            "Chromosome-held-out\n(LOCO, honest)")),
  auroc = c(prim$random_split_auroc, prim$auroc)
)
gap <- prim$random_split_auroc - prim$auroc  # 0.0591

p2 <- ggplot(d2, aes(x = split, y = auroc, group = 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3,
             color = "#9E9E9E") +
  geom_line(linewidth = 0.4, color = "#9E9E9E") +
  geom_point(aes(color = split), size = 2.4) +
  geom_text(aes(label = sprintf("%.3f", auroc)),
            vjust = -1.1, size = GEOM_TEXT_6PT, color = "black") +
  annotate("text", x = 1.5, y = mean(d2$auroc) + 0.004,
           label = sprintf("gap %.3f\nLD leakage, not biology", gap),
           size = GEOM_TEXT_6PT, color = "black", lineheight = 0.9) +
  scale_color_manual(values = c("Random split\n(leaky)" = "#9E9E9E",
                                "Chromosome-held-out\n(LOCO, honest)" = masld_colors$up),
                     guide = "none") +
  scale_y_continuous(limits = c(0.49, 0.60), breaks = seq(0.50, 0.60, 0.02)) +
  labs(x = NULL, y = "Sign-concordance auROC",
       title = "Cross-validation leakage gap") +
  theme_masld()

message("[leakage_gap] CAPTION: Supervised direction model evaluated under ",
        "random-split CV (0.567, leaky) vs chromosome-held-out LOCO (0.507, honest). ",
        "The 0.059 optimism is LD leakage between train/test SNPs in the same locus, ",
        "not a real biological signal; dashed line = 0.5 chance.")

save_fig(p2, file.path(OUT_DIR, "leakage_gap.pdf"),
         width = 2.8, height = 2.4)
rendered <- c(rendered, "leakage_gap.pdf")

# =============================================================================
# PANEL 3 — coding_gof_lof.pdf
# =============================================================================
# 42 coding effectors as a per-model call heatmap. Columns: AlphaMissense,
# ESM1b, popEVE (missense severity) + LoGoFunc (GoF/LoF/Neutral) + Curated
# literature. Spotlight the 3 neomorphs (PNPLA3 I148M / TM6SF2 E167K / APOE
# C130R) that every sequence/function predictor calls benign/Neutral -> the
# curated GoF flags must be retained. Held-out GoF-vs-LoF auROC 0.905 annotated.
# -----------------------------------------------------------------------------
cod <- fread(file.path(SF, "coding_hardening_v2.tsv"))

# Harmonize each model's native call to a shared discrete scale.
sev_map <- c(
  AM_lbenign = "Benign/Tolerant", AM_ambiguous = "Moderate/Ambiguous",
  AM_lpath   = "Damaging/Pathogenic",
  ESM_tolerant = "Benign/Tolerant", ESM_moderate = "Moderate/Ambiguous",
  ESM_damaging = "Damaging/Pathogenic",
  popEVE_tolerant = "Benign/Tolerant", popEVE_moderate = "Moderate/Ambiguous",
  popEVE_damaging = "Damaging/Pathogenic"
)
lgf_map <- c(Neutral = "Neutral", LOF = "LoF", GOF = "GoF", GoF = "GoF")
# Curated-literature flag (parsed from gof_lof_flag / gof_lof_call).
# Curated column shares the function-class labels (no duplicate legend swatches);
# the neomorph distinction is carried by the marker + caption, not a legend entry.
cod[, curated := fifelse(grepl("GoF", gof_lof_call, ignore.case = TRUE) &
                           grepl("neomorph", gof_lof_call, ignore.case = TRUE),
                         "GoF",
                 fifelse(gof_lof_call == "LoF", "LoF",
                 fifelse(gof_lof_call == "isoform_functional", "Isoform",
                         "none")))]

cod[, `:=`(
  AlphaMissense = sev_map[am_call],
  ESM1b         = sev_map[esm_call],
  popEVE        = sev_map[popeve_call],
  LoGoFunc      = lgf_map[logofunc_call]
)]
cod[is.na(AlphaMissense), AlphaMissense := NA_character_]

# Spotlight the 3 canonical lipid-genetics effectors at the top; only PNPLA3 I148M
# is a true GoF neomorph. TM6SF2 E167K = curated LoF, APOE C130R = isoform-functional.
spotlight <- c("PNPLA3 I148M", "TM6SF2 E167K", "APOE C130R")
neomorphs <- c("PNPLA3 I148M")
cod[, glab := paste(gene, protein_variant)]
cod[, is_spot := glab %in% spotlight]
cod[, is_neo  := glab %in% neomorphs]

# Density trim: drop the all-neutral/all-benign rows that carry no signal (they made the
# grid mostly-empty and buried the message). Keep any variant with a non-benign severity
# call, a LoF/GoF discriminator call, a curated flag, or a spotlight effector.
sev_hit <- c("Moderate/Ambiguous", "Damaging/Pathogenic")
cod <- cod[is_spot | curated != "none" | LoGoFunc %in% c("LoF", "GoF") |
             AlphaMissense %in% sev_hit | ESM1b %in% sev_hit | popEVE %in% sev_hit]

# Long form. Column order left->right; keep curated as the right-most strip.
mcols <- c("AlphaMissense", "ESM1b", "popEVE", "LoGoFunc", "Curated")
long <- rbindlist(list(
  cod[, .(glab, gene, is_neo, model = "AlphaMissense", call = AlphaMissense)],
  cod[, .(glab, gene, is_neo, model = "ESM1b",         call = ESM1b)],
  cod[, .(glab, gene, is_neo, model = "popEVE",        call = popEVE)],
  cod[, .(glab, gene, is_neo, model = "LoGoFunc",      call = LoGoFunc)],
  cod[, .(glab, gene, is_neo, model = "Curated",       call = curated)]
))
long[call %in% c("", "none"), call := NA_character_]
long[, model := factor(model, levels = mcols)]

# Row order: LoGoFunc LoF block first, then Neutral; spotlight effectors pulled to top.
ord <- cod[order(-is_spot, factor(logofunc_call, levels = c("LOF","Neutral","")),
                 gene, protein_variant), glab]
long[, glab := factor(glab, levels = rev(ord))]

# One clean discrete scale: missense severity (3) + function class (LoF/GoF/Isoform/Neutral).
# No duplicate swatches; curated cells reuse the same function-class colors.
fill_vals <- c(
  "Benign/Tolerant"     = "#E8E8E8",  # light gray
  "Moderate/Ambiguous"  = "#F4A674",  # Liang warm peach
  "Damaging/Pathogenic" = "#C9265E",  # Liang deep magenta
  "Neutral"             = "#9E9E9E",  # control gray
  "LoF"                 = "#1565C0",  # deep blue
  "GoF"                 = "#7B1FA2",  # violet
  "Isoform"             = "#00695C"   # teal
)
long[, call := factor(call, levels = names(fill_vals))]

# Italic gene symbol + plain variant for y labels; dagger marks neomorphs.
ylabs <- setNames(levels(long$glab), levels(long$glab))
neo_set <- ord[ord %in% neomorphs]

p3 <- ggplot(long, aes(x = model, y = glab, fill = call)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 4.5, linewidth = 0.35, color = "black") +  # sep predictors|curated
  scale_fill_manual(values = fill_vals, na.value = "#FFFFFF",
                    name = NULL, drop = TRUE) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Coding effectors: per-model function calls") +
  theme_masld() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
        axis.text.y = element_text(size = 5, face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.22, "cm"))

# Mark the 3 neomorph rows with a magenta point at the left margin + annotate auROC.
neo_df <- data.frame(glab = factor(neo_set, levels = rev(ord)),
                     model = factor("AlphaMissense", levels = mcols))
p3 <- p3 +
  annotate("point", x = 0.35,
           y = which(rev(ord) %in% neo_set), size = 1.7, shape = 18,
           color = "#7B1FA2") +
  coord_cartesian(xlim = c(0.3, 5.5), clip = "off")

message("[coding_gof_lof] CAPTION: 42 curated coding effectors, per-model ",
        "function calls. Missense-severity predictors (AlphaMissense, ESM1b, ",
        "popEVE) and the LoGoFunc GoF/LoF/Neutral discriminator; right strip = ",
        "curated-literature flag. The violet dot (left) marks the sole GoF neomorph, ",
        "PNPLA3 I148M -- called Neutral by LoGoFunc (P_GoF 0.042); no predictor emits ",
        "a GoF call, so the curated GoF flag is retained. TM6SF2 E167K (curated LoF) ",
        "and APOE C130R (curated isoform-functional) are the other two spotlighted ",
        "lipid effectors, curated to their true mechanisms. LoGoFunc held-out ",
        "GoF-vs-LoF auROC = 0.905 [0.879,0.930], perm-p 5e-4 (variant-level held-out ",
        "upper bound; not gene-family disjoint).")

save_fig(p3, file.path(OUT_DIR, "coding_gof_lof.pdf"),
         width = 3.6, height = 6.6)
rendered <- c(rendered, "coding_gof_lof.pdf")

# =============================================================================
# PANEL 4 — qtl_modality_direction.pdf   [LANDED 2026-07-11]
# =============================================================================
# The KEY new panel: direction is learnable by MODALITY. The same supervised
# elastic-net (nested LOCO, pooled-OOF) is at chance for eQTL EXPRESSION direction
# (0.507) but predictive for caQTL ACCESSIBILITY direction (0.707, beats both
# baselines CI-separated, perm-p 1e-3, leakage-gap ~0). sQTL (GTEx) is
# leakage-flagged (in model training) -> uncertified, shown gray with asterisk.
# GUARDRAIL: bounded benchmark; caQTL/sQTL never a scored convergence channel;
# NOT a paired test (different variant sets & n) -- framed as "accessibility
# direction is learnable where expression direction is not".
# -----------------------------------------------------------------------------
vj  <- jsonlite::fromJSON(file.path(SF, "direction_hardened_verdict.json"))
pan <- vj$panels
mkrow <- function(nm, disp, verdict) data.table(
  panel   = disp,
  auroc   = pan[[nm]]$enet_auroc,
  ci_lo   = pan[[nm]]$enet_ci_lo,
  ci_hi   = pan[[nm]]$enet_ci_hi,
  n       = pan[[nm]]$n,
  tier    = pan[[nm]]$leakage_tier,
  verdict = verdict)
d4 <- rbindlist(list(
  mkrow("caqtl", "caQTL — accessibility", "Predictable (beats baselines)"),
  mkrow("eqtl",  "eQTL — expression",     "At chance"),
  mkrow("sqtl",  "sQTL — splice*",        "Leakage-flagged (uncertified)")
))
d4[, panel := factor(panel, levels = d4[order(auroc), panel])]

p4 <- ggplot(d4, aes(x = auroc, y = panel, color = verdict)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3,
             color = "#9E9E9E") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.16, linewidth = 0.35) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c(
      "Predictable (beats baselines)" = masld_colors$up,
      "At chance"                     = "#9E9E9E",
      "Leakage-flagged (uncertified)" = "#9E9E9E"),
    name = NULL) +
  scale_x_continuous(limits = c(0.35, 0.76), breaks = seq(0.4, 0.75, 0.1)) +
  labs(x = "Sign-concordance auROC (pooled OOF)", y = NULL,
       title = "Direction is learnable by modality") +
  theme_masld() +
  theme(legend.position = "bottom", legend.margin = margin(t = -2))

message("[qtl_modality_direction] CAPTION: Same supervised elastic-net ",
        "(nested leave-one-chromosome-out, pooled-OOF sign-auROC) applied to three ",
        "QTL modalities. caQTL ACCESSIBILITY direction is predictable (0.707 ",
        "[0.697,0.718], perm-p 1e-3, beats aggregated-signed and positional baselines ",
        "CI-separated, leakage-gap ~0, n=11,896), whereas eQTL EXPRESSION direction ",
        "sits at chance (0.507 [0.431,0.577], n=211). sQTL (GTEx, in model training) is ",
        "leakage-flagged and uncertified (asterisk). Point estimates on different ",
        "variant sets -- NOT a paired test. Interpretation: sequence-to-function models resolve the ",
        "proximal sequence->chromatin layer but not emergent expression direction. ",
        "Bounded supplementary benchmark; APPLY-ONLY -- never a scored convergence channel.")

save_fig(p4, file.path(OUT_DIR, "qtl_modality_direction.pdf"),
         width = 3.6, height = 2.2)
rendered <- c(rendered, "qtl_modality_direction.pdf")

# =============================================================================
# PANEL 5 — decima_celltype.pdf   [LANDED 2026-07-11]
# =============================================================================
# Cell-type expression VEP (Decima): argmax liver compartment of the 120
# nomination/coding/lead variants. ~1/3 hepatocyte, ~2/3 non-parenchymal ->
# ties the genetic arm to the multi-cellular cascade. Hepatocyte highlighted;
# non-parenchymal in control gray. GUARDRAIL: magnitude/mechanism annotation
# only, never a scored channel; snATAC corroboration is semi-independent.
# -----------------------------------------------------------------------------
dec <- fread(file.path(SF, "decima_celltype.tsv"))
d5  <- dec[, .N, by = argmax_celltype][order(-N)]
d5[, pretty := gsub("_", " ", argmax_celltype)]
d5[, is_hep := argmax_celltype == "hepatocyte"]
d5[, pretty := factor(pretty, levels = rev(pretty))]
n_tot <- d5[, sum(N)]
pct_hep <- round(100 * d5[is_hep == TRUE, sum(N)] / n_tot)

p5 <- ggplot(d5, aes(x = N, y = pretty, fill = is_hep)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = N), hjust = -0.3, size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(values = c("TRUE" = masld_colors$up, "FALSE" = "#9E9E9E"),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Variants (argmax cell type)", y = NULL,
       title = "Cell-type localization of genetic nominations") +
  theme_masld()

message("[decima_celltype] CAPTION: Decima cell-type EXPRESSION VEP argmax ",
        sprintf("compartment for %d nomination/coding/lead variants: %d%% hepatocyte, ",
                n_tot, pct_hep),
        sprintf("%d%% non-parenchymal (macrophage/lymphoid/LSEC/stellate). ", 100 - pct_hep),
        "SORT1 rs12740374 -> hepatocyte (positive control PASS). snATAC peak-presence ",
        "agrees 19/29 (semi-independent, coverage-limited -- same Human-Multiome cohort, ",
        "NOT fully disjoint). Hepatocyte highlighted. Magnitude/mechanism annotation ",
        "only; corroborates the multi-cellular cascade; never a scored convergence channel.")

save_fig(p5, file.path(OUT_DIR, "decima_celltype.pdf"),
         width = 3.4, height = 1.9)
rendered <- c(rendered, "decima_celltype.pdf")

# =============================================================================
cat("\n=== Phase-5 supplementary figures rendered ===\n")
for (f in rendered) {
  fp <- file.path(OUT_DIR, f)
  cat(sprintf("  %-26s %s  (%d bytes)\n", f,
              ifelse(file.exists(fp), "OK", "MISSING"),
              ifelse(file.exists(fp), file.info(fp)$size, 0L)))
}
cat("Output dir:", OUT_DIR, "\n")
