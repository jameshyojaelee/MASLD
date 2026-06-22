#!/usr/bin/env Rscript
# =============================================================================
# fig_disease_signatures_presentation.R
# Presentation-ready overview of all bulk RNA-seq disease signature analyses.
# 8 standalone panels for individual slides.
#
# Output: figures/presentation/disease_signature/panel_{A-H}.pdf
# Env: micromamba activate rnaseq
# SLURM: --partition=cpu --mem=16G --cpus-per-task=4 --time=48:00:00
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Presentation constants ───────────────────────────────────────────────────
PB   <- 13        # base font
TXT  <- 4.2       # primary label
TXTs <- 3.6       # secondary label
TXTt <- 3.0       # tertiary label
PT   <- 3.0       # point size
LW   <- 0.7       # primary linewidth
LWt  <- 0.5       # thin linewidth

OUT <- file.path(BASE, "figures/presentation/disease_signature")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

spf <- function(p, fname, w, h) {
  dir.create(dirname(fname), recursive = TRUE, showWarnings = FALSE)
  dev <- if (capabilities("cairo")) cairo_pdf else pdf
  ggsave(fname, p, width = w, height = h, device = dev)
}

cat("=== Disease Signature Presentation Figures ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# Path shortcuts
INT_RES  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
PROG_DIR <- file.path(INT_RES, "progression")
MPROG    <- file.path(INT_RES, "multiprogram")
SUB_DIR  <- file.path(BASE, "RNA-seq/results/subtypes")

# =============================================================================
# Panel A: Dream Mega-Analysis Volcano
# =============================================================================
cat("--- Panel A: volcano ---\n")

dream_path <- file.path(INT_RES, "integration/canonical_deg_results.csv")
dream <- fread(dream_path)

# canonical_deg_results.csv uses unprefixed logFC/padj (C2 cutover 2026-06-08)
if ("adj.P.Val" %in% names(dream) && !"padj" %in% names(dream))
  setnames(dream, "adj.P.Val", "padj")
stopifnot(all(c("padj", "logFC") %in% names(dream)))

dream[, sig := fifelse(padj < 0.05 & abs(logFC) > 0.5,
  fifelse(logFC > 0, "Up", "Down"), "NS")]
dream[, neg_log_p := pmin(-log10(padj), 50)]  # cap for display

n_up   <- sum(dream$sig == "Up")
n_down <- sum(dream$sig == "Down")
n_deg  <- n_up + n_down

# Top genes for labeling
dream[, abs_t := abs(get(names(dream)[grep("^t$|^t\\.", names(dream))[1]]))]
top_genes <- dream[sig != "NS"][order(-abs_t)][1:min(15, .N)]

pA <- ggplot(dream, aes(x = logFC, y = neg_log_p, color = sig)) +
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 0.4), dpi = 300) +
  geom_point(data = top_genes, size = 1.5, alpha = 0.9) +
  geom_text_repel(data = top_genes, aes(label = symbol),
    size = TXTt, max.overlaps = 20, segment.size = 0.3,
    min.segment.length = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
    linewidth = LWt, color = "gray60") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed",
    linewidth = LWt, color = "gray60") +
  scale_color_manual(values = c(Up = masld_colors$up,
    Down = masld_colors$down, NS = masld_colors$ns), guide = "none") +
  annotate("text", x = max(dream$logFC) * 0.7, y = 48,
    label = sprintf("%s DEGs\n(%s up, %s down)",
      comma(n_deg), comma(n_up), comma(n_down)),
    size = TXTs, color = "gray25", fontface = "italic") +
  labs(x = expression(log[2]*" fold change"),
       y = expression(-log[10]*"(padj)"),
       title = "Dream mega-analysis (5 control-bearing cohorts, 847 samples)") +
  theme_masld(base_size = PB) +
  theme(plot.title = element_text(lineheight = 1.1))

spf(pA, file.path(OUT, "panel_A_volcano.pdf"), w = 9, h = 6)
cat("  Written: panel A\n")

# =============================================================================
# Panel B: Progression Cascade (DEGs per fibrosis transition)
# =============================================================================
cat("--- Panel B: progression cascade ---\n")

trans_path <- file.path(PROG_DIR, "transition_summary.csv")
if (file.exists(trans_path)) {
  trans <- fread(trans_path)
  trans_fib <- trans[stage_type == "fibrosis"]

  # Build stacked up/down
  stack_dt <- rbindlist(list(
    trans_fib[, .(transition, direction = "Upregulated",
                  count = n_up)],
    trans_fib[, .(transition, direction = "Downregulated",
                  count = -n_down)]
  ))
  # Clean labels
  stack_dt[, transition := gsub("_to_", " \u2192 ", transition)]
  fib_order <- c("F0 \u2192 F1", "F1 \u2192 F2", "F2 \u2192 F3", "F3 \u2192 F4")
  stack_dt[, transition := factor(transition, levels = fib_order)]

  pB <- ggplot(stack_dt, aes(x = transition, y = count, fill = direction)) +
    geom_col(width = 0.65) +
    geom_hline(yintercept = 0, linewidth = LWt) +
    geom_text(aes(label = comma(abs(count)),
                  y = count + sign(count) * 150),
      size = TXTs, color = "gray25") +
    scale_fill_manual(values = c(Upregulated = masld_colors$up,
      Downregulated = masld_colors$down), name = NULL) +
    labs(x = "Fibrosis stage transition",
         y = "Number of DEGs (padj < 0.05)",
         title = "Progression cascade across fibrosis transitions") +
    theme_masld(base_size = PB) +
    theme(legend.position = "top")

  spf(pB, file.path(OUT, "panel_B_progression_cascade.pdf"), w = 8, h = 5.5)
  cat("  Written: panel B\n")
} else {
  cat("  SKIP: transition_summary.csv not found\n")
}

# =============================================================================
# Panel C: Mid-stage (F1-F3) Inflection Classification
# =============================================================================
cat("--- Panel C: mid-stage (F1-F3) inflection classification ---\n")

switch_path <- file.path(BASE, "RNA-seq/results/stratified_causal/switch_gene_classification.csv")
if (file.exists(switch_path)) {
  sw <- fread(switch_path)

  # Count by classification
  sw_counts <- sw[, .N, by = classification]
  sw_counts[, group := fifelse(grepl("switch", classification), "Switch-like", "Gradual")]
  sw_counts[, direction := fifelse(grepl("up", classification), "Up", "Down")]

  group_tot <- sw_counts[, .(total = sum(N), pct = sum(N) / nrow(sw) * 100), by = group]

  # Bar chart
  sw_counts[, label := paste0(gsub("_", " ", classification), "\n(n=", comma(N), ")")]
  sw_counts[, classification := factor(classification,
    levels = c("switch_up", "switch_down", "gradual_up", "gradual_down"))]

  class_colors <- c(
    switch_up   = "#E91E63",  # bright magenta
    switch_down = "#42A5F5",  # light blue
    gradual_up  = "#C2185B",  # deep magenta
    gradual_down = "#1565C0"  # deep blue
  )
  class_labels <- c(
    switch_up   = "Switch \u2191",
    switch_down = "Switch \u2193",
    gradual_up  = "Gradual \u2191",
    gradual_down = "Gradual \u2193"
  )

  sw_counts[, display := class_labels[as.character(classification)]]
  sw_counts[, display := factor(display, levels = class_labels)]

  n_switch <- sw_counts[group == "Switch-like", sum(N)]
  pct_switch <- round(n_switch / nrow(sw) * 100, 1)

  pC <- ggplot(sw_counts, aes(x = display, y = N, fill = classification)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = comma(N), y = N + 60), size = TXTs) +
    scale_fill_manual(values = class_colors, guide = "none") +
    annotate("text", x = 1.5, y = max(sw_counts$N) * 0.9,
      label = sprintf("%s switch-like genes\n(%.1f%% of DEGs)", comma(n_switch), pct_switch),
      size = TXT, color = "gray25", fontface = "bold") +
    labs(x = NULL, y = "Number of DEGs",
         title = "Mid-stage (F1-F3) inflection classification (AIC: step vs linear model)") +
    theme_masld(base_size = PB)

  spf(pC, file.path(OUT, "panel_C_midstage_inflection.pdf"), w = 8, h = 5.5)
  cat("  Written: panel C\n")
} else {
  cat("  SKIP: switch_gene_classification.csv not found\n")
}

# =============================================================================
# Panel D: Consensus Tiers
# =============================================================================
cat("--- Panel D: consensus tiers ---\n")

atlas_path <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
if (file.exists(atlas_path)) {
  atlas <- fread(atlas_path, select = c("human_consensus_tier"))
  tier_dt <- atlas[!is.na(human_consensus_tier) & human_consensus_tier != "",
    .N, by = human_consensus_tier]
  tier_dt[, tier := human_consensus_tier]

  tier_dt[, tier := factor(tier, levels = c("Tier1_HighConfidence",
    "Tier2_Moderate", "Tier3_Exploratory"))]
  tier_dt[, label := gsub("_", " ", as.character(tier))]
  tier_dt[, label := factor(label, levels = gsub("_", " ", levels(tier_dt$tier)))]

  pD <- ggplot(tier_dt, aes(x = label, y = N, fill = tier)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = comma(N), y = N + 10), size = TXT) +
    scale_fill_manual(values = tier_consensus_colors, guide = "none") +
    labs(x = NULL, y = "Number of genes",
         title = "Consensus DEG tiers (dream + per-cohort replication)") +
    theme_masld(base_size = PB) +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))

  spf(pD, file.path(OUT, "panel_D_consensus_tiers.pdf"), w = 7, h = 5)
  cat("  Written: panel D\n")
} else {
  cat("  SKIP: multi_evidence_atlas.csv not found\n")
}

# =============================================================================
# Panel E: Conserved (Cross-Species)
# =============================================================================
cat("--- Panel E: conserved core ---\n")

conc_path <- file.path(BASE,
  "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv")
if (file.exists(conc_path)) {
  conc <- fread(conc_path, select = c("primary_category"))

  # Map to display categories
  cat_map <- c(
    Conserved       = "Conserved",
    Human_Enriched       = "Human-specific",
    Species_Discordant   = "Species-discordant",
    Mouse_Specific       = "Mouse-specific",
    Moderate_Concordance = "Other",
    Diet_Selective       = "Other",
    Not_Significant      = "Other"
  )
  conc[, display_cat := cat_map[primary_category]]
  conc[is.na(display_cat), display_cat := "Other"]

  cat_counts <- conc[, .N, by = display_cat]
  cat_counts[, display_cat := factor(display_cat,
    levels = c("Conserved", "Human-specific", "Mouse-specific",
               "Species-discordant", "Other"))]

  disp_colors <- c(
    "Conserved"     = concordance_colors[["Conserved"]],
    "Human-specific"     = concordance_colors[["Human-specific"]],
    "Mouse-specific"     = concordance_colors[["Mouse-specific"]],
    "Species-discordant" = concordance_colors[["Species-discordant"]],
    "Other"              = concordance_colors[["Other"]]
  )

  pE <- ggplot(cat_counts, aes(x = display_cat, y = N, fill = display_cat)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = comma(N), y = N + 150), size = TXTs) +
    scale_fill_manual(values = disp_colors, guide = "none") +
    labs(x = NULL, y = "Number of genes",
         title = "Cross-species concordance (human vs 5 mouse diet models)") +
    theme_masld(base_size = PB) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  spf(pE, file.path(OUT, "panel_E_conserved.pdf"), w = 9, h = 5.5)
  cat("  Written: panel E\n")
} else {
  cat("  SKIP: concordance_atlas_unified.csv not found\n")
}

# =============================================================================
# Panel F: Deconvolution Attribution
# =============================================================================
cat("--- Panel F: deconvolution attribution ---\n")

deconv_path <- file.path(BASE, "RNA-seq/results/causal_inference/deconv_attribution_scores.csv")
if (file.exists(deconv_path)) {
  deconv <- fread(deconv_path)

  # Only significant genes
  deconv_sig <- deconv[category != "Not_significant"]
  attr_counts <- deconv_sig[, .N, by = category]
  attr_counts[, pct := round(N / sum(N) * 100, 1)]
  attr_counts[, label := paste0(gsub("_", " ", category), "\n(", comma(N), ", ", pct, "%)")]
  attr_counts[, category := factor(category,
    levels = c("Hepatocyte_intrinsic", "Composition_driven", "Unmasked"))]

  pF <- ggplot(attr_counts, aes(x = category, y = N, fill = category)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = paste0(comma(N), "\n(", pct, "%)"),
                  y = N + 100), size = TXTs, lineheight = 0.9) +
    scale_fill_manual(values = attribution_colors, guide = "none") +
    scale_x_discrete(labels = c("Hepatocyte\nintrinsic",
                                "Composition\ndriven", "Unmasked")) +
    labs(x = NULL, y = "Number of DEGs",
         title = "Deconvolution attribution of disease DEGs") +
    theme_masld(base_size = PB)

  spf(pF, file.path(OUT, "panel_F_attribution.pdf"), w = 7, h = 5)
  cat("  Written: panel F\n")
} else {
  cat("  SKIP: deconv_attribution_scores.csv not found\n")
}

# =============================================================================
# Panel G: Multi-Program Ordinal Decomposition
# =============================================================================
cat("--- Panel G: multiprogram ablation ---\n")

abl_path <- file.path(MPROG, "multiprogram_ablation.csv")
if (file.exists(abl_path)) {
  abl <- fread(abl_path)

  # Clean task and modality names
  task_map <- c(fibrosis = "Fibrosis\n(F0\u2013F4)",
                nas_composite = "NAS\n(0\u20138)",
                severity = "Severity\n(4-class)")
  mod_map <- c(expression = "Expression", celltype = "Cell-type",
               ssgsea = "ssGSEA", coloc_genetics = "COLOC",
               tf_activity = "TF activity", pseudotime = "Pseudotime",
               clinical = "Clinical", combined = "Combined")

  abl[, task_label := task_map[task]]
  abl[, mod_label  := mod_map[modality]]
  abl <- abl[!is.na(task_label) & !is.na(mod_label)]

  abl[, task_label := factor(task_label, levels = unname(task_map))]
  abl[, mod_label  := factor(mod_label,  levels = unname(mod_map))]

  pG <- ggplot(abl, aes(x = mod_label, y = mean_qwk, fill = task_label)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, mean_qwk - std_qwk),
                      ymax = pmin(1, mean_qwk + std_qwk)),
      position = position_dodge(width = 0.75), width = 0.2,
      linewidth = LWt, color = "gray40") +
    scale_fill_manual(values = c("#C2185B", "#7B1FA2", "#1565C0"),
      name = "Target") +
    labs(x = NULL, y = "Quadratic Weighted Kappa (QWK)",
         title = "Multi-program ordinal prediction by feature group") +
    theme_masld(base_size = PB) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "top")

  spf(pG, file.path(OUT, "panel_G_multiprogram.pdf"), w = 10, h = 5.5)
  cat("  Written: panel G\n")
} else {
  cat("  SKIP: multiprogram_ablation.csv not found\n")
}

# =============================================================================
# Panel H: NMF Subtyping (de-emphasized)
# =============================================================================
cat("--- Panel H: NMF subtyping ---\n")

nmf_assign  <- fread(file.path(SUB_DIR, "nmf_assignments.csv"))
nmf_metrics <- file.path(SUB_DIR, "nmf_metrics.csv")
clin_path   <- file.path(SUB_DIR, "subtype_clinical_associations.csv")
meta_path   <- file.path(INT_RES, "../../metadata/unified_metadata.csv")

# k-program refactor (2026-04-20): nmf_subtype is a legacy binary alias where
# S2 = dominant program is Fibrotic, S1 = any non-Fibrotic program lumped.
# Relabel for presentation.
nmf_assign[, subtype_label := fifelse(nmf_subtype == "S1", "Non-Fibrogenic", "Fibrogenic")]

if (file.exists(meta_path)) {
  meta <- fread(meta_path, select = c("sample_id", "fibrosis_stage"))
  comp <- merge(nmf_assign, meta, by = "sample_id")
  comp <- comp[!is.na(fibrosis_stage)]
  comp[, fib_label := paste0("F", fibrosis_stage)]

  prop_dt <- comp[, .N, by = .(fib_label, subtype_label)]
  prop_dt[, total := sum(N), by = fib_label]
  prop_dt[, pct := N / total * 100]
  prop_dt[, fib_label := factor(fib_label, levels = paste0("F", 0:4))]

  fib_pct <- prop_dt[subtype_label == "Fibrogenic"]
  setorder(fib_pct, fib_label)

  pH <- ggplot(prop_dt, aes(x = fib_label, y = pct, fill = subtype_label)) +
    geom_col(width = 0.65) +
    geom_text(data = fib_pct,
              aes(y = pct / 2, label = sprintf("%.0f%%", pct)),
              color = "white", size = TXTs, fontface = "bold") +
    scale_fill_manual(values = c(Metabolic = masld_colors$down,
                                 Fibrogenic = masld_colors$up),
      name = "NMF subtype") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "Fibrosis stage", y = "Subtype proportion",
         title = "NMF molecular subtypes across fibrosis stages") +
    theme_masld(base_size = PB) +
    theme(legend.position = "top")

  spf(pH, file.path(OUT, "panel_H_nmf_subtyping.pdf"), w = 7, h = 5)
  cat("  Written: panel H\n")
} else {
  cat("  SKIP: unified_metadata.csv not found\n")
}

cat("\n=== Disease Signature Presentation Figures Complete ===\n")
cat("End:", format(Sys.time()), "\n")
cat("Outputs in:", OUT, "\n")
