#!/usr/bin/env Rscript
# =============================================================================
# Figure 1: Atlas Overview (revised per FIGURE_PLAN_REVISED.md, 2026-04-15)
#
# Panels:
#   1a: Study-design schematic — supplied externally, NOT generated here.
#   1b: Sunburst (fig_sunburst.py) — supplied externally, NOT generated here.
#   1c: BulkRNA metadata matrix (fig_bulkrna_matrix.py --panel) — NOT generated here.
#   1d: UpSet: per-study DEGs vs integrated DEGs (padj<0.05 primary ~12,266; ashr lfsr<0.05 supp. ~6,501; STAR -s 2 2026-05-28; counts computed at runtime)
#   1e: Per-study LOO-CV replication (rho = 0.959, recovery 88.1% [62-95%]) —
#       rendered to figS01/panels/figS01_loo_recovery.pdf (demoted to supp)
#   1f: Per-study replication of integrated DEGs (ggplot RDS produced by
#       fig_cohort_replication.R; promoted from figS01)
#   1g: LOO-CV stability of integrated DEGs vs |log2FC| cutoff (heatmap +
#       across-fold mean-recovery curve; ggplot RDS produced by
#       fig1g_loo_lfc_stability.R; motivates the LFC cutoff via stabilization)
#
# Outputs:
#   figures/main/fig1_atlas_overview/panels/fig1{d,f,g}.pdf
#   figures/supplementary/figS_methods_validation/qc_validation/panels/figS01_loo_recovery.pdf
#   (Combined fig1_atlas_overview.pdf is NOT generated — panels are arranged
#    manually.)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG1_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Cohort name map (accession) ----
STUDY_NAMES <- c(
  GSE126848   = "GSE126848",
  GSE130970   = "GSE130970",
  GSE135251   = "GSE135251",
  GSE162694   = "GSE162694",
  GSE167523   = "GSE167523",
  GSE174478   = "GSE174478",
  GSE193066   = "GSE193066",
  GSE213621   = "GSE213621",
  GSE240729   = "GSE240729"
)
# PRJNA512027 (Gerhard 2018) excluded from cohort presentation: L0/S0
# library-prep batch perfectly confounded with diagnosis (all 34 controls
# L0, all 102 NASH S0). Still loaded by the pipeline for the
# fibrosis-vs-healthy contrast.

# Stage palette — fig1_colors (palette1/2/3-derived; softer than masld_colors)
stage_palette <- c(
  Control       = fig1_colors$control,
  NAFL          = fig1_colors$nafl,
  Borderline    = fig1_colors$borderline,
  NASH          = fig1_colors$nash,
  Fibrosis_only = fig1_colors$fibrosis,
  Unknown       = fig1_colors$unknown
)

# =============================================================================
# Load data once
# =============================================================================
cat("Loading QC + metadata...\n")
qc_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")
qc <- fread(qc_path)
meta_path <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
meta <- fread(meta_path)
meta <- merge(meta, qc[, .(sample_id, pass_technical)], by = "sample_id", all.x = TRUE)
meta_qc <- meta[pass_technical == TRUE]
# Drop PRJNA512027 from cohort presentation (L0/S0 batch confound — see
# top-of-file note). 1,444 raw QC-pass -> 1,259 after PRJNA512027 removal.
meta_qc <- meta_qc[dataset != "PRJNA512027"]
cat(sprintf("  QC-pass samples (9 cohorts, PRJNA excluded): %d\n", nrow(meta_qc)))

# Stage harmonization
meta_qc[, stage := diagnosis_harmonized]
meta_qc[is.na(stage) | stage == "", stage := NA_character_]
# Fibrosis-only cohorts (no NASH/NAFL labels) use fibrosis_stage if available
fibrosis_only_datasets <- c("GSE213621", "GSE240729")
meta_qc[dataset %in% fibrosis_only_datasets & is.na(stage) & !is.na(fibrosis_stage),
        stage := "Fibrosis_only"]
meta_qc[is.na(stage), stage := "Unknown"]
# Group Borderline with NASH per NAS>=3 harmonization (Kleiner 2005)
# Keep Borderline distinct for color but annotate
meta_qc[, stage := factor(stage, levels = c("Control", "NAFL", "Borderline",
                                             "NASH", "Fibrosis_only", "Unknown"))]

# Cohort totals
cohort_tot <- meta_qc[, .(n = .N), by = dataset][order(-n)]
cohort_tot[, author := STUDY_NAMES[dataset]]

# Panel 1a is the study-design schematic supplied externally (not generated here).

# PANEL 1b: Cohort x stage stacked bar (1,259 QC-pass denominator, 9 cohorts)
#   Horizontal stacks: one bar per cohort, stacked by diagnosis_harmonized.
#   Width proportional to cohort QC-pass size (area still maps to sample count).
# =============================================================================
cat("\n── Panel 1b: Cohort × stage stacked bar ──\n")

tm <- meta_qc[, .N, by = .(dataset, stage)]
tm[, author := STUDY_NAMES[dataset]]
tm[, cohort_total := sum(N), by = author]
tm[, frac := N / cohort_total]
tm[, author_label := sprintf("%s\n(n=%d)", author, cohort_total)]
# Sort cohorts by total size (largest on top)
cohort_order <- tm[, .(tot = cohort_total[1]), by = author_label][order(tot)]
tm[, author_label := factor(author_label, levels = cohort_order$author_label)]

# Stage stack order
tm[, stage := factor(stage, levels = c("Control", "NAFL", "Borderline",
                                        "NASH", "Fibrosis_only", "Unknown"))]

# Use absolute counts so bar length encodes cohort N.
p_b <- ggplot(tm, aes(x = N, y = author_label, fill = stage)) +
  geom_col(color = "white", linewidth = 0.3, width = 0.8) +
  geom_text(data = tm[N >= 15],
            aes(label = N),
            position = position_stack(vjust = 0.5),
            color = "white", size = PUB_GEOM_TEXT, fontface = "plain") +
  scale_fill_manual(values = stage_palette,
                    name = "NAS-harmonised\ndiagnosis",
                    drop = FALSE,
                    labels = c(Control = "Control",
                               NAFL = "NAFL (MASL)",
                               Borderline = "Borderline (NAS 3\u20134)",
                               NASH = "NASH (MASH)",
                               Fibrosis_only = "Fibrosis-only label",
                               Unknown = "Unlabelled")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)),
                     labels = comma,
                     breaks = pretty_breaks(n = 5)) +
  labs(x = "QC-pass samples", y = NULL) +
  theme_masld() + theme_pub() +
  theme(legend.position = "right")
message(sprintf("[caption] 9 cohorts \u00b7 %s QC-pass samples \u00b7 diagnosis composition",
                format(nrow(meta_qc), big.mark = ",")))

# DO NOT WRITE fig1b.pdf — the canonical 1b is the convergence-wheel sunburst
# (scripts/figures/fig_sunburst.py). Cohort × stage stack moved to figS01.
# save_fig(p_b, file.path(PANEL_DIR, "fig1b.pdf"),
#          width = fig_full_width, height = 3.4)
# cat("  Saved fig1b.pdf\n")

# =============================================================================
# PANEL 1c: Metadata availability grid per cohort
#  Rows: cohorts; cols: sample count, fibrosis stage, NAS, sex, age
# =============================================================================
cat("\n── Panel 1c: Metadata overview grid ──\n")

avail <- meta_qc[, .(
  n_samples          = .N,
  n_fibrosis         = sum(!is.na(fibrosis_stage) & fibrosis_stage != ""),
  n_nas              = sum(!is.na(nas_score) & nas_score != ""),
  n_sex              = sum(!is.na(sex) & sex != ""),
  n_age              = sum(!is.na(age) & age != "")
), by = dataset]
avail[, author := STUDY_NAMES[dataset]]
avail <- avail[order(-n_samples)]
avail[, author := factor(author, levels = rev(author))]

# Fractions
avail[, frac_fib := n_fibrosis / n_samples]
avail[, frac_nas := n_nas / n_samples]
avail[, frac_sex := n_sex / n_samples]
avail[, frac_age := n_age / n_samples]

long <- melt(avail,
             id.vars = c("author", "n_samples"),
             measure.vars = c("frac_fib", "frac_nas", "frac_sex", "frac_age"),
             variable.name = "field", value.name = "frac")
long_counts <- melt(avail,
                    id.vars = c("author", "n_samples"),
                    measure.vars = c("n_fibrosis", "n_nas", "n_sex", "n_age"),
                    variable.name = "field", value.name = "count")
long_counts[, field := factor(field,
                              levels = c("n_fibrosis", "n_nas", "n_sex", "n_age"),
                              labels = c("frac_fib", "frac_nas", "frac_sex", "frac_age"))]
long <- merge(long, long_counts[, .(author, field, count)],
              by = c("author", "field"))

field_labels <- c(frac_fib = "Fibrosis\nstage",
                  frac_nas = "NAS\nscore",
                  frac_sex = "Sex",
                  frac_age = "Age")
long[, field := factor(field,
                       levels = c("frac_fib", "frac_nas", "frac_sex", "frac_age"),
                       labels = field_labels)]

long[, label := ifelse(frac > 0,
                        sprintf("%d/%d", count, n_samples),
                        "\u2014")]

# Plot: heatmap (availability) + left-side sample-count bar
p_c_heat <- ggplot(long, aes(x = field, y = author, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label),
            color = ifelse(long$frac > 0.55, "white", "black"),
            size = PUB_GEOM_TEXT) +
  scale_fill_gradient(low = "#F5F5F5", high = fig1_colors$up,
                      limits = c(0, 1), labels = percent_format(accuracy = 1),
                      name = "Annotation\navailable") +
  labs(x = NULL, y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(angle = 0))
message("[caption] Per-cohort metadata availability — 615/892 fibrosis-annotated samples used for C1 fibrosis-adjusted contrast")

p_c_bar <- ggplot(avail, aes(x = n_samples, y = author)) +
  geom_col(fill = fig1_colors$down, width = 0.7) +
  geom_text(aes(label = format(n_samples, big.mark = ",")),
            hjust = -0.15, size = PUB_GEOM_TEXT, color = "gray20") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)),
                     labels = comma) +
  labs(x = "QC-pass samples", y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p_c <- p_c_heat + p_c_bar + plot_layout(widths = c(2.2, 1.1))

# DO NOT WRITE fig1c.pdf — the canonical 1c is the bulkRNA metadata matrix
# (scripts/figures/fig_bulkrna_matrix.py --panel). Availability heatmap moved
# to figS01.
# save_fig(p_c, file.path(PANEL_DIR, "fig1c.pdf"),
#          width = fig_full_width, height = 3.0)
# cat("  Saved fig1c.pdf\n")

# =============================================================================
# PANEL 1d: UpSet — per-study DEGs vs integrated DEGs
#  Headline: primary DEGs (padj<0.05, |logFC|>0.5); ashr supplementary
# =============================================================================
cat("\n── Panel 1d: UpSet per-study vs integrated DEGs ──\n")

dream <- load_dream_results()
per_study <- load_per_study_de()

# Use disease-vs-control contrast studies (fibrosis-only studies use rank contrast)
upset_studies <- c("GSE126848", "GSE130970", "GSE135251",
                   "GSE162694", "GSE213621")
upset_names <- STUDY_NAMES[upset_studies]

# Threshold: primary (padj<0.05, |logFC|>0.5)
padj_thr <- 0.05
lfc_thr  <- 0.5

deg_lists <- list()
for (ds in upset_studies) {
  ds_data <- per_study[dataset == ds]
  pcol <- intersect(c("padj", "adj.P.Val"), names(ds_data))[1]
  lcol <- intersect(c("logFC", "study_logFC"), names(ds_data))[1]
  if (is.na(pcol) || is.na(lcol)) next
  degs <- ds_data[get(pcol) < padj_thr & abs(get(lcol)) > lfc_thr, gene]
  degs <- sub("\\..*", "", degs)
  deg_lists[[upset_names[ds]]] <- unique(degs)
}
dream_degs <- dream[is_dream_deg(dream), gene]
dream_degs <- sub("\\..*", "", dream_degs)
deg_lists[["Integrated"]] <- unique(dream_degs)

set_names <- names(deg_lists)
all_genes_upset <- unique(unlist(deg_lists))
bin_dt <- data.table(gene = all_genes_upset)
for (s in set_names) bin_dt[, (s) := as.integer(gene %in% deg_lists[[s]])]
bin_dt[, int_key := do.call(paste0, .SD), .SDcols = set_names]
int_counts <- bin_dt[, .N, by = int_key][order(-N)]
# Curated intersection set: top 10 by size + the "Integrated only" intersection
# pinned in (it sits well below the top-10 by count, but is the headline panel
# story — genes recovered only by the meta-analysis).
int_only_key <- paste0(strrep("0", length(set_names) - 1), "1")
int_counts[, n_active := nchar(gsub("0", "", int_key))]
# Curated set: top 10 by size + "Integrated only" (always pinned) +
# every high-replication intersection (active in >=5 of the 6 sets), so the
# 5-/6-way bars are visible even when their counts fall below the top-10 cut.
top_ints <- int_counts[1:min(10, nrow(int_counts))]
extra_keys <- int_counts[(n_active >= 5 & N > 3) | int_key == int_only_key, int_key]
extras <- int_counts[int_key %in% extra_keys & !int_key %in% top_ints$int_key]
if (nrow(extras) > 0) top_ints <- rbind(top_ints, extras)
top_ints[, n_active := NULL]
top_ints <- top_ints[order(-N)]
top_n <- nrow(top_ints)
top_ints[, rank := .I]
for (i in seq_along(set_names)) {
  top_ints[, (set_names[i]) := as.integer(substr(int_key, i, i))]
}
set_sizes <- vapply(set_names, function(s) sum(bin_dt[[s]]), integer(1))
set_labels_vec <- paste0(set_names, " (", format(set_sizes, big.mark = ","), ")")

dot_long <- melt(top_ints, id.vars = c("int_key", "N", "rank"),
                 measure.vars = set_names,
                 variable.name = "set", value.name = "active")
dot_long[, set := factor(set, levels = set_names)]
seg_data <- dot_long[active == 1, .(ymin = min(as.numeric(set)),
                                    ymax = max(as.numeric(set))), by = rank]
seg_data <- seg_data[ymin != ymax]

top_ints[, has_dream := as.logical(get("Integrated"))]
n_active <- rowSums(as.matrix(top_ints[, ..set_names]))
# Three distinct fills: integration-only (magenta), mixed integrated+per-study
# (violet), per-study only (blue). Earlier scheme reused masld_colors$up for
# both integration-only and mixed bars, making them visually identical.
top_ints[, bar_fill := fifelse(has_dream & n_active == 1, fig1_colors$up,
                         fifelse(has_dream, fig1_colors$mixed, fig1_colors$down))]

p_d_bars <- ggplot(top_ints, aes(x = rank, y = N)) +
  geom_col(fill = top_ints$bar_fill, width = 0.85) +
  geom_text(aes(label = format(N, big.mark = ",")),
            vjust = -0.3, size = PUB_GEOM_TEXT) +
  scale_x_continuous(limits = c(0.4, top_n + 0.6), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)), labels = comma) +
  labs(y = "Intersection\nsize") +
  theme_masld() + theme_pub() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = margin(2, 5, 0, 5))
message(sprintf("[caption] Per-study vs integrated DEGs \u00b7 %s integrated DEGs (effect-size-aware interval-null FDR gate, FDR<0.05 at lfc=0.25)",
                format(set_sizes["Integrated"], big.mark = ",")))

p_d_dots <- ggplot() +
  geom_segment(data = seg_data,
               aes(x = rank, xend = rank, y = ymin, yend = ymax),
               color = "gray30", linewidth = 0.4) +
  geom_point(data = dot_long,
             aes(x = rank, y = as.numeric(set), fill = factor(active)),
             shape = 21, size = 2, stroke = 0.3, color = "gray40") +
  scale_fill_manual(values = c("0" = "#E8E8E8", "1" = "gray20"),
                    guide = "none") +
  scale_x_continuous(limits = c(0.4, top_n + 0.6), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(set_names),
                     labels = set_labels_vec,
                     limits = c(0.5, length(set_names) + 0.5),
                     expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0, 5, 2, 5))

p_d <- p_d_bars / p_d_dots + plot_layout(heights = c(2, 1.4))

# fig1d.pdf relocated to the Figure-3 RNA-seq dir (FIG2_DIR = figures/main/fig3_RNAseq,
# back-compat constant name). All OTHER panels in this script still write to FIG1_DIR.
save_fig(p_d, file.path(FIG2_DIR, "panels", "cohort_diagnosis_composition.pdf"),
         width = fig_col_width, height = 2.6)
cat("  Saved cohort_diagnosis_composition.pdf\n")

# =============================================================================
# PANEL 1e: Per-study LOO-CV replication
# =============================================================================
cat("\n── Panel 1e: LOO-CV per-study replication ──\n")

LOO_DIR <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/loo_cv")
loo <- fread(file.path(LOO_DIR, "loo_cv_summary.csv"))
loo[, author := STUDY_NAMES[held_out]]
# Cohort N (QC-pass)
loo <- merge(loo,
             qc[pass_technical == TRUE, .(n_samples = .N), by = dataset],
             by.x = "held_out", by.y = "dataset", all.x = TRUE)
loo[, author := factor(author, levels = author[order(pct_full_recovered)])]

mean_rec <- mean(loo$pct_full_recovered)
rng_rec <- range(loo$pct_full_recovered)
mean_rho <- mean(loo$spearman_rho)
worst <- loo[which.min(pct_full_recovered)]

p_e_rec <- ggplot(loo, aes(x = pct_full_recovered, y = author)) +
  geom_segment(aes(x = 0, xend = pct_full_recovered,
                   y = author, yend = author),
               color = fig1_colors$ns, linewidth = 0.4) +
  geom_point(aes(size = n_samples), color = fig1_colors$up, shape = 16) +
  geom_vline(xintercept = mean_rec, linetype = "dashed",
             color = "gray30", linewidth = 0.4) +
  annotate("text", x = mean_rec, y = 0.6,
           label = sprintf("mean %.1f%%", mean_rec),
           size = PUB_GEOM_TEXT, color = "gray20", hjust = -0.1) +
  scale_size_continuous(range = c(2, 5), name = "N") +
  scale_x_continuous(limits = c(55, 105),
                     breaks = c(60, 70, 80, 90, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "DEG recovery (% of full model)",
       y = "Held-out cohort") +
  theme_masld() + theme_pub() +
  theme(legend.position = c(0.9, 0.2))
message(sprintf("[caption] Per-study LOO-CV \u00b7 range %.0f-%.0f%%", rng_rec[1], rng_rec[2]))

p_e_rho <- ggplot(loo, aes(x = spearman_rho, y = author)) +
  geom_segment(aes(x = 0.75, xend = spearman_rho,
                   y = author, yend = author),
               color = fig1_colors$ns, linewidth = 0.4) +
  geom_point(aes(size = n_samples), color = fig1_colors$down, shape = 16) +
  geom_vline(xintercept = mean_rho, linetype = "dashed",
             color = "gray30", linewidth = 0.4) +
  annotate("text", x = mean_rho, y = 0.6,
           label = sprintf("mean \u03c1 = %.3f", mean_rho),
           size = PUB_GEOM_TEXT, color = "gray20", hjust = -0.05) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  scale_x_continuous(limits = c(0.75, 1.01),
                     breaks = c(0.8, 0.9, 1.0)) +
  labs(x = "Spearman \u03c1 vs full model",
       y = NULL) +
  theme_masld() + theme_pub() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
message(sprintf("[caption] GSE213621 weakest fold \u00b7 %d%% / rho = %.2f",
                round(worst$pct_full_recovered), worst$spearman_rho))

p_e <- p_e_rec + p_e_rho + plot_layout(widths = c(1.6, 1))

# Panel E (per-study LOO-CV replication) lives in figS01 — it was removed
# from the main fig1 layout, so write it directly to the supplementary panels
# directory rather than the main-fig panels dir.
SUPP_PANEL_DIR <- file.path(FIGS01_DIR, "panels")
dir.create(SUPP_PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
save_fig(p_e, file.path(SUPP_PANEL_DIR, "figS01_loo_recovery.pdf"),
         width = fig_full_width, height = 3.0)
cat("  Saved figS01/panels/figS01_loo_recovery.pdf\n")

# =============================================================================
# PANEL 1f: Integrated-DEG cross-cohort replication (promoted from figS01)
# =============================================================================
# Built by scripts/figures/fig_cohort_replication.R, which serializes the
# ggplot to FIGS01_DIR/cohort_replication_pC.rds and writes fig1f.pdf to
# this PANEL_DIR. We re-load the ggplot here so it composes into the combined
# patchwork rather than just sitting alongside as an orphan PDF.
cat("\n── Panel 1f: integrated-DEG cross-cohort replication ──\n")

pC_rds <- file.path(FIGS01_DIR, "cohort_replication_pC.rds")
if (file.exists(pC_rds)) {
  p_f <- readRDS(pC_rds)
  has_p_f <- TRUE
  cat("  Loaded p_f from ", pC_rds, "\n", sep = "")
} else {
  warning("fig1f source RDS not found: ", pC_rds,
          " -- run fig_cohort_replication.R first. Combined fig1 will skip panel F.")
  has_p_f <- FALSE
}

# =============================================================================
# PANEL 1g: Patient-consistency vs |log2FC| heatmap, Up+Down combined
# =============================================================================
cat("\n── Panel 1g: patient-consistency × LFC cutoff heatmap ──\n")

pG_rds <- file.path(FIGS01_DIR, "fig1g_loo_stability_pG.rds")
if (file.exists(pG_rds)) {
  p_g <- readRDS(pG_rds)
  has_p_g <- TRUE
  cat("  Loaded p_g from ", pG_rds, "\n", sep = "")
} else {
  warning("fig1g source RDS not found: ", pG_rds,
          " -- run fig1g_loo_lfc_stability.R first. Combined fig1 will skip panel G.")
  has_p_g <- FALSE
}

# Combined fig1_atlas_overview.pdf is intentionally NOT generated. Panels are
# assembled by hand from figures/main/fig1_atlas_overview/panels/*.pdf.

cat("\n=== Figure 1 complete (sub-panels only) ===\n")
