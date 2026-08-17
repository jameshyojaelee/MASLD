#!/usr/bin/env Rscript
# KEY MESSAGE: Adjacent-stage expression burden and continuous bulk program activity describe cross-sectional remodeling along fibrosis and NAS axes.
# ============================================================================
# progression_cascade.R  — Figure 3E (cross-sectional stage remodeling)
#
# Two-column layout sharing three matched track types:
#   LEFT  column — fibrosis stage (Kleiner F0..F4), built inline here
#   RIGHT column — NAS activity (NAS0..NAS5-8), reused from nas_activity_cascade.R
# Track types (rows), row labels on the fibrosis (left) column only:
#   TRACK 1 — stage-specific DEG counts (up vs down) vs the axis baseline
#   TRACK 2 — NMF 6-program dominant-sample composition shift
#   TRACK 3 — cell-type composition shift (hepatocyte decline, NPC rise)
#
# The NAS tracks are sourced (read-only) from nas_activity_cascade.R rather than
# re-derived, so the two columns stay definition-identical. Restores the pre-
# 2026-07-02 combined NAS|fibrosis cascade as the main panel (fibrosis-left).
#
# Output: FIG2_DIR/panels/fig3e_progression_cascade.pdf  (Fig 3E; combined axes)
#   (FIG2_DIR resolves to figures/main/fig3_bulk_transcriptomics.)
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# Candidate six-figure layout. This branch uses the synchronized fragment-count
# release and continuous prespecified P1-P6 axes. The legacy renderer below is
# retained only for reproducibility of the existing main-figure PDF.
CANDIDATE_ROOT <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
if (nzchar(CANDIDATE_ROOT)) {
  STAGE_ROOT <- normalizePath(Sys.getenv("STAGE_RELEASE_ROOT", ""), mustWork = TRUE)
  EXT_ROOT <- file.path(CANDIDATE_ROOT, "analysis", "stage_extensions")
  OUT_DIR <- file.path(CANDIDATE_ROOT, "figure3", "panels")
  SUPP_DIR <- file.path(CANDIDATE_ROOT, "supplementary", "figureS3", "panels")
  SRC_DIR <- file.path(CANDIDATE_ROOT, "source_tables")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)

  adjacent <- fread(file.path(EXT_ROOT, "stage_all_gene_results.tsv"))
  extension <- fread(file.path(EXT_ROOT, "stage_extension_all_gene_results.tsv"))
  activity <- fread(file.path(EXT_ROOT, "prespecified_nmf_program_activity.tsv"))
  audit <- fread(file.path(EXT_ROOT, "stage_extension_design_audit.tsv"))
  expected_family <- unique(audit$bh_family_size)
  if (length(expected_family) != 1L ||
      any(adjacent[, uniqueN(gene_id_versioned), by = contrast]$V1 != expected_family) ||
      any(audit$bh_family_size != expected_family)) stop("Incomplete stage BH family", call. = FALSE)

  # The main fibrosis bars use the same F0 reference as the centered NMF
  # activity tracks. Adjacent-stage contrasts remain available in the source
  # release but are not mixed with the common-reference display.
  fibrosis_vs_f0 <- extension[axis == "fibrosis"]
  fibrosis_vs_f0[, display := comparison]
  nas <- extension[axis == "NAS"]
  nas[, display := comparison]
  counts <- rbind(
    fibrosis_vs_f0[FDR < 0.05,
      .(n = .N), by = .(axis, display, direction = fifelse(logFC > 0, "Higher", "Lower"))],
    data.table(axis = "fibrosis", display = "F0",
               direction = c("Higher", "Lower"), n = 0L),
    nas[FDR < 0.05,
      .(n = .N), by = .(axis, display, direction = fifelse(logFC > 0, "Higher", "Lower"))],
    data.table(axis = "NAS", display = "NAS0",
               direction = c("Higher", "Lower"), n = 0L)
  )
  counts[, signed_n := fifelse(direction == "Lower", -n, n)]
  counts[, direction := factor(direction, levels = c("Higher", "Lower"))]

  count_plot <- function(axis_id, groups, xlab, ylab, labels = groups) {
    d <- copy(counts[axis == axis_id])
    d[, display := factor(display, levels = groups, labels = labels)]
    # A group whose contrast returns a handful of genes draws a bar too short to
    # see against a several-thousand-gene axis, which reads as missing data
    # rather than as a near-null step. Print the count for any non-zero bar
    # under 2% of the panel maximum. The definitional zero at the reference
    # group stays unlabelled.
    tiny <- d[n > 0 & n < 0.02 * max(abs(d$signed_n))]
    ggplot(d, aes(display, signed_n, fill = direction)) +
      geom_col(width = 0.44) +
      geom_hline(yintercept = 0, linewidth = 0.25) +
      geom_text(data = tiny, aes(label = n,
                                 vjust = ifelse(direction == "Higher", -0.5, 1.4)),
                size = 6 / .pt, show.legend = FALSE) +
      scale_fill_manual(values = c(Higher = masld_colors$mash, Lower = masld_colors$down), guide = "none") +
      scale_y_continuous(labels = abs, expand = expansion(mult = c(0.13, 0.16))) +
      labs(x = xlab, y = ylab) + theme_masld_compact() +
      theme(axis.text.x = element_text(size = 6, face = "plain"),
            plot.margin = margin(2, 2, 1, 2))
  }
  program_colors <- setNames(cat_palette[seq_len(6)], paste0("P", 1:6))
  program_plot <- function(axis_id, groups, labels, xlab) {
    d <- activity[axis == axis_id]
    d[, group := factor(group, levels = groups, labels = labels)]
    ggplot(d, aes(group, centered_activity, group = program, color = program)) +
      geom_hline(yintercept = 0, color = "grey75", linewidth = 0.25) +
      geom_line(linewidth = 0.55) + geom_point(size = 1.05) +
      scale_color_manual(values = program_colors, guide = guide_legend(nrow = 2)) +
      # drop = FALSE keeps a group with no estimate as an empty tick so the two
      # rows of a column stay on the same x geometry.
      scale_x_discrete(drop = FALSE) +
      labs(x = xlab, y = "NMF activity", color = NULL) + theme_masld_compact() +
      theme(axis.text.x = element_text(size = 6, face = "plain"),
            legend.position = "bottom", legend.key.width = unit(8, "pt"),
            plot.margin = margin(1, 2, 2, 2))
  }
  fib_count <- count_plot("fibrosis", paste0("F", 0:4), "Fibrosis stage", "DEGs vs F0")
  # NAS uses the same common-reference estimand as the fibrosis column: each
  # activity group against strict NAS0, with NAS0 itself the definitional zero.
  # The adjacent NAS chain stays on axis "NAS_adjacent" for the supplement.
  nas_count <- count_plot("NAS", c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8"),
                          "NAS group", "DEGs vs NAS0",
                          labels = c("0", "1–2", "3–4", "5–8"))
  fib_program <- program_plot("fibrosis", paste0("F", 0:4), paste0("F", 0:4), "Fibrosis stage")
  # Only one NAS0 donor on the synchronized fragment substrate has a
  # prespecified NMF score, so program activity begins at NAS1–2 and is centred
  # there, not on NAS0. The NAS0 tick is carried as an empty position to hold
  # the column geometry; do not imply a supported NAS0 program estimate. The
  # DEG track above is unaffected and uses all 46 eligible NAS0 donors.
  nas_program <- program_plot("NAS", c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8"),
                              c("0", "1–2", "3–4", "5–8"), "NAS group")
  # Main panel is two matched columns: fibrosis (F0 + 4 groups) left, NAS
  # (NAS0 + 3 groups) right. Flat 2x2 so plot_layout() widths apply to the real
  # columns; the NAS column is proportionally narrower for its one fewer group
  # rather than stretched to match. The right program plot drops its duplicated
  # y title only (ticks stay) because both program tracks are in the same
  # centered NMF-activity units.
  main <- fib_count + nas_count +
    fib_program + (nas_program + labs(y = NULL)) +
    plot_layout(ncol = 2, heights = c(1, 1.05), widths = c(1, 0.85),
                guides = "collect") &
    theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "pt"))
  ggsave(file.path(OUT_DIR, "fig3e_stage_remodeling.pdf"), main,
         width = 5.00, height = 2.27, device = cairo_pdf)

  # Composition remains supplementary. Each biopsy is one participant-level
  # unit under the validated crosswalk; no library or cell is treated as n.
  comp <- fread(file.path(
    BASE,
    "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/",
    "workstreams/BULK-PROGRAM-MAP-v9/precontract/composition_acceptance/accepted_composition.tsv.gz"
  ))
  stage_samples <- unique(fread(file.path(EXT_ROOT, "stage_sample_manifest.tsv"))[, .(
    sample_id, analysis_unit_id, fibrosis_stage
  )], by = "analysis_unit_id")
  stage_samples[, group := paste0("F", fibrosis_stage)]
  nas_samples <- fread(file.path(EXT_ROOT, "stage_extension_sample_manifest.tsv"))[axis == "NAS"]
  nas_samples[, group := fifelse(grepl("NAS0", contrast) & group == "NAS0", "NAS0", group)]
  composition_summary <- function(samples, axis_id, groups) {
    d <- merge(comp, unique(samples[, .(sample_id, analysis_unit_id, group)], by = "analysis_unit_id"), by = "sample_id")
    cell_cols <- setdiff(names(comp), c("sample_id", "dataset"))
    d[, hepatocyte := Hepatocytes]
    d[, non_parenchymal := rowSums(.SD, na.rm = TRUE), .SDcols = setdiff(cell_cols, "Hepatocytes")]
    z <- melt(d, id.vars = c("analysis_unit_id", "group"), measure.vars = c("hepatocyte", "non_parenchymal"),
              variable.name = "compartment", value.name = "fraction")
    z[, axis := axis_id]
    z[, .(mean = mean(fraction), SE = sd(fraction) / sqrt(.N), n_participants = .N),
      by = .(axis, group = factor(group, levels = groups), compartment)]
  }
  comp_sum <- rbind(
    composition_summary(stage_samples, "fibrosis", paste0("F", 0:4)),
    composition_summary(nas_samples, "NAS", c("NAS0", "NAS1-2", "NAS3-4", "NAS5-8"))
  )
  comp_plot <- function(axis_id, xlab) ggplot(comp_sum[axis == axis_id],
      aes(group, mean, color = compartment, group = compartment)) +
    geom_line(linewidth = 0.65) + geom_point(size = 1.1) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = 0.1, linewidth = 0.3) +
    scale_color_manual(values = c(hepatocyte = "#9E9E9E", non_parenchymal = cat_palette[[3]]),
                       labels = c("Hepatocyte", "Non-parenchymal"), name = NULL) +
    labs(x = xlab, y = "Estimated fraction") + theme_masld_compact() +
    theme(axis.text.x = element_text(size = 6, face = "plain"), legend.position = "bottom")
  full <- (fib_count | nas_count) / (fib_program | nas_program) /
    (comp_plot("fibrosis", "Fibrosis stage") | comp_plot("NAS", "NAS group")) +
    plot_layout(heights = c(1, 1.05, 0.9), guides = "collect") & theme(legend.position = "bottom")
  ggsave(file.path(SUPP_DIR, "figs3_stage_remodeling_full.pdf"), full,
         width = 7.1, height = 4.8, device = cairo_pdf)

  fwrite(counts, file.path(SRC_DIR, "fig3e_deg_counts.tsv"), sep = "\t")
  fwrite(activity, file.path(SRC_DIR, "fig3e_prespecified_nmf_activity.tsv"), sep = "\t")
  fwrite(comp_sum, file.path(SRC_DIR, "figs3_composition_summary.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}

PANEL_DIR <- file.path(FIG2_DIR, "panels")
DATA_DIR  <- file.path(PANEL_DIR, "data")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR,  showWarnings = FALSE, recursive = TRUE)

STAGES <- c("F0", "F1", "F2", "F3", "F4")
PANEL_WIDTH_IN  <- 4.40
PANEL_HEIGHT_IN <- 2.27

# DEG-calling mode for Track 1. CANONICAL = treat @ lfc=0.25 — the repo-wide
# TREAT definition (treat_fdr<0.05 @ lfc=0.25; McCarthy & Smyth 2009), which
# folds the effect-size floor INTO the significance test. This matches the NAS
# cascade and the rest of the paper. Env overrides exist for sensitivity only:
#   "treat"   — TREAT interval-null FDR at lfc=LFC (default, canonical).
#   "sig"     — significance only: padj < FDR, NO effect-size floor.
#   "default" — CANONICAL from 2026-08-12: padj < FDR & |raw log2FC| > LFC.
# Canonical gate migrated 2026-08-12 from TREAT (lfc=0.25) to the conventional
# adjusted-p + effect-size floor. The floor applies to the RAW logFC, the same
# statistic the paper-wide canonical uses; shrunk_logFC is kept for direction only.
DEG_METHOD <- tolower(Sys.getenv("PROG_CASCADE_DEG_METHOD", "default"))
FDR_CUTOFF <- as.numeric(Sys.getenv("PROG_CASCADE_FDR", "0.05"))
LFC_CUTOFF <- as.numeric(Sys.getenv("PROG_CASCADE_LFC", "0.5"))
stopifnot(DEG_METHOD %in% c("treat", "sig", "default"))
# The canonical run writes the plain back-compat filename; sensitivity arms get a tag.
is_canonical <- DEG_METHOD == "default" && LFC_CUTOFF == 0.5 && FDR_CUTOFF == 0.05
file_tag <- if (is_canonical) "" else switch(DEG_METHOD,
  sig     = "_sig",
  treat   = sprintf("_treat_lfc%s", gsub("[.]", "p", format(LFC_CUTOFF, trim = TRUE))),
  default = sprintf("_legacy_lfc%s", gsub("[.]", "p", format(LFC_CUTOFF, trim = TRUE))))

# Shared x-axis geometry and margins keep F0..F4 columns aligned vertically
# across all tracks. Cell-type labels live in the collected legend rather than
# extending the bottom panel beyond its x-scale.
X_EXPAND      <- ggplot2::expansion(mult = c(0.04, 0.04))
SHARED_MARGIN <- ggplot2::margin(3, 3, 3, 3)

# ----------------------------------------------------------------------------
# TRACK 1 — stage-specific DEG counts (up vs down) vs F0
# ----------------------------------------------------------------------------
deg <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_dream.csv"))
# Effect size: prefer ashr-shrunk; fall back to raw logFC if shrunk missing.
deg[, eff := if ("shrunk_logFC" %in% names(deg)) shrunk_logFC else logFC]
deg[is.na(eff), eff := logFC]
deg[, stage := sub("_vs_F0$", "", contrast)]              # F1_vs_F0 -> F1
if (DEG_METHOD == "treat") {
  # Reconstruct limma::treat from the saved moderated coef/SE/t/P value; df.total
  # is common within each fitted contrast (mirrors nas_activity_cascade.R).
  infer_df <- function(t_stat, p_value) {
    ok <- which(is.finite(t_stat) & t_stat != 0 & is.finite(p_value) &
                p_value > 0 & p_value < 1)[1]
    uniroot(function(df) 2 * pt(-abs(t_stat[ok]), df) - p_value[ok],
            c(1, 1e7), tol = 1e-10)$root
  }
  deg[, treat_df := infer_df(t, P.Value), by = contrast]
  deg[, treat_p := { a <- abs(logFC)
    pt((a - LFC_CUTOFF) / SE, df = treat_df, lower.tail = FALSE) +
      pt((a + LFC_CUTOFF) / SE, df = treat_df, lower.tail = FALSE) }]
  deg[, treat_fdr := p.adjust(treat_p, method = "BH"), by = contrast]
  deg_sig <- deg[treat_fdr < FDR_CUTOFF]
} else if (DEG_METHOD == "sig") {
  deg_sig <- deg[padj < FDR_CUTOFF]                        # significance only
} else {
  # CANONICAL 2026-08-12: floor on the RAW logFC, matching is_canonical_deg().
  deg_sig <- deg[padj < FDR_CUTOFF & abs(logFC) > LFC_CUTOFF]
}
deg_counts <- deg_sig[, .(
    up   = sum(eff > 0, na.rm = TRUE),
    down = sum(eff < 0, na.rm = TRUE)
  ), by = stage]
# F0 vs F0 = baseline, zero by definition (no contrast computed)
deg_counts <- rbind(
  data.table(stage = "F0", up = 0L, down = 0L),
  deg_counts, fill = TRUE)
deg_counts[, stage := factor(stage, levels = STAGES)]
setorder(deg_counts, stage)

deg_long <- melt(deg_counts, id.vars = "stage",
                 measure.vars = c("up", "down"),
                 variable.name = "direction", value.name = "n")
# Down plotted as negative for a back-to-back (diverging) read.
deg_long[, n_signed := ifelse(direction == "down", -n, n)]
deg_long[, direction := factor(direction, levels = c("up", "down"),
                               labels = c("Up vs F0", "Down vs F0"))]

DEG_COLORS <- c("Up vs F0"   = masld_colors$mash,   # #C9265E magenta
                "Down vs F0" = masld_colors$down)   # #1565C0 blue

p1 <- ggplot(deg_long, aes(stage, n_signed, fill = direction)) +
  geom_col(width = 0.66, color = NA) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  geom_text(data = deg_long[n > 0],
            aes(label = n,
                vjust = ifelse(direction == "Up vs F0", -0.35, 1.25)),
            size = GEOM_TEXT_6PT, color = "black") +
  scale_fill_manual(values = DEG_COLORS, name = NULL) +
  scale_y_continuous(labels = function(x) abs(x), expand = expansion(mult = c(0.12, 0.18))) +
  scale_x_discrete(expand = X_EXPAND) +
  labs(x = NULL, y = "Stage DEGs\n(vs F0)", fill = "Stage DEG") +
  theme_masld_compact() +
  theme(legend.key.size = unit(0.28, "cm"),
        legend.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = SHARED_MARGIN)

# ----------------------------------------------------------------------------
# TRACK 2 — NMF program composition: fraction of samples dominated by each P1..P6
# ----------------------------------------------------------------------------
nmf  <- fread(file.path(BASE, "RNA-seq/results/subtypes/nmf_assignments.csv"))
labs <- fread(file.path(BASE, "RNA-seq/results/subtypes/program_labels.csv"))
meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"),
  select = c("sample_id", "fibrosis_stage"))

# program_labels uses program_code/biological_label (not program/label).
lab_map <- setNames(labs$biological_label, labs$program_code)

nmf <- merge(nmf[, .(sample_id, dominant_program_code)],
             meta, by = "sample_id")
# Exclude coarse-staged GSE213621 (F0F1/F3F4 grouped) + dropped PRJNA512027 so the
# F0-F4 axis is true Kleiner, matching Track 1's DEG source. See load_figure_data.R.
clean_staged_ids <- staged_disease_meta(disease_only = FALSE)$sample_id
nmf <- nmf[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4 & sample_id %in% clean_staged_ids]
nmf[, stage := factor(paste0("F", fibrosis_stage), levels = STAGES)]
nmf[, prog_label := lab_map[dominant_program_code]]

nmf_comp <- nmf[, .N, by = .(stage, dominant_program_code, prog_label)]
nmf_comp[, frac := N / sum(N), by = stage]
# Order programs from disease/active (top of stack) to quiescent/stable (bottom).
prog_order <- c("P1", "P3", "P2", "P4", "P5", "P6")
prog_labels_ordered <- lab_map[prog_order]
nmf_comp[, prog_label := factor(prog_label, levels = prog_labels_ordered)]

# Colors: active programs get saturated cat hues; the two quiescent-parenchyma
# baselines and the two stable programs anchor in/near control gray.
NMF_COLORS <- c(
  "Inflammatory-EMT"     = cat_palette[3],  # red-coral — active disease      # was: "Progression-Inflammatory"
  "Fibrotic-ECM"         = cat_palette[6],  # violet    — active fibrosis      # was: "Fibrogenic"
  "Hepatocyte-Metabolic" = "#BDBDBD",       # near-gray (baseline parenchyma)  # was: "Quiescent-Parenchyma_1"
  "Unresolved"           = "#D6D6D6",       #                                   # was: "Quiescent-Parenchyma_2"
  "Noncoding"            = "#9E9E9E",       # control gray                      # was: "Stable_1"
  "Skeletal-muscle"      = "#7D7D7D")       #                                   # was: "Stable_2"

p2 <- ggplot(nmf_comp, aes(stage, frac, fill = prog_label)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = NMF_COLORS, name = "NMF program",
                    breaks = prog_labels_ordered,
                    labels = prog_labels_ordered) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = ""),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(expand = X_EXPAND) +
  labs(x = NULL, y = "NMF program\n(% samples)") +
  # 2026-07-09 resize (4.40x2.27in target): ncol=2 trades legend WIDTH for
  # HEIGHT (6 items: 3 rows+title instead of 6+title) -- at fixed 6pt font the
  # single-column legend stack no longer fits the shorter target height.
  guides(fill = guide_legend(ncol = 2, keyheight = unit(0.22, "cm"))) +
  theme_masld_compact() +
  theme(legend.position = "right",
        legend.text  = element_text(size = 6, face = "plain"),
        legend.title = element_text(size = 6, face = "plain"),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        plot.margin  = SHARED_MARGIN)

# ----------------------------------------------------------------------------
# TRACK 3 — cell-type composition shift across stages
# ----------------------------------------------------------------------------
ct <- fread(file.path(BASE,
  "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv"))
ct <- ct[!is.na(fibrosis_stage) & fibrosis_stage %in% 0:4 & sample_id %in% clean_staged_ids]
ct[, stage := factor(paste0("F", fibrosis_stage), levels = STAGES)]

CT_KEEP <- c("Hepatocytes", "Fibroblasts", "Macrophages",
             "T cells", "Endothelial cells")
ct_mean <- ct[, lapply(.SD, mean, na.rm = TRUE), by = stage,
              .SDcols = CT_KEEP]
ct_long <- melt(ct_mean, id.vars = "stage",
                variable.name = "celltype", value.name = "frac")
# Sort cell types by direction of change F0->F4: hepatocyte (declining) first,
# then the rising non-parenchymal compartments by magnitude of rise.
delta <- ct_long[, .(d = frac[stage == "F4"] - frac[stage == "F0"]), by = celltype]
ct_order <- c("Hepatocytes",
              as.character(delta[celltype != "Hepatocytes"][order(-d)]$celltype))
ct_long[, celltype := factor(celltype, levels = ct_order)]

# Hepatocyte baseline = control gray (it is the compartment being displaced);
# rising NPC compartments take semantic ct_palette hues.
CT_COLORS <- c(
  "Hepatocytes"       = "#9E9E9E",        # control gray (displaced baseline)
  "Fibroblasts"       = ct_palette[["Fibroblasts"]],
  "Macrophages"       = ct_palette[["Macrophages"]],
  "T cells"           = ct_palette[["T cells"]],
  "Endothelial cells" = ct_palette[["Endothelial cells"]])

# Hepatocyte (78-93%) and the NPC compartments (0-4%) live on incompatible
# scales; a single linear axis crushes the NPC rise. Split Track 3 into two
# sub-panels that still share the F0..F4 x-axis: (3a) hepatocyte decline,
# (3b) the rising non-parenchymal compartments (magnified, NOT % of total ->
# % within compartment is misleading; keep absolute fraction with its own axis).
ct_hep <- ct_long[celltype == "Hepatocytes"]
ct_npc <- ct_long[celltype != "Hepatocytes"]
ct_npc[, celltype := droplevels(celltype)]

p3a <- ggplot(ct_hep, aes(stage, frac, group = 1)) +
  geom_line(linewidth = 0.8, color = "#9E9E9E") +
  geom_point(size = 1.4, color = "#9E9E9E") +
  geom_text(data = ct_hep[stage %in% c("F0", "F4")],
            aes(label = scales::percent(frac, accuracy = 1),
                hjust = ifelse(stage == "F0", -0.2, 1.2)),
            vjust = -0.5, size = GEOM_TEXT_6PT, color = "black") +
  coord_cartesian(clip = "off") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = ""),
                     limits = c(min(ct_hep$frac) * 0.97, max(ct_hep$frac) * 1.06)) +
  scale_x_discrete(expand = X_EXPAND) +
  labs(x = NULL, y = "Hepatocyte\n(%)") +
  theme_masld_compact() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = SHARED_MARGIN)

p3b <- ggplot(ct_npc, aes(stage, frac, color = celltype, group = celltype)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  scale_color_manual(values = CT_COLORS, name = "Cell type",
                     breaks = levels(ct_npc$celltype)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = ""),
                     expand = expansion(mult = c(0.05, 0.10))) +
  scale_x_discrete(expand = X_EXPAND) +
  labs(x = "Fibrosis stage (Kleiner F0–F4)", y = "Non-parenchymal\n(%)") +
  # 2026-07-09 resize: ncol=2 (4 items: 2 rows+title instead of 4+title).
  guides(color = guide_legend(ncol = 2, keyheight = unit(0.22, "cm"))) +
  theme_masld_compact() +
  theme(legend.position = "right",
        legend.text = element_text(size = 6, face = "plain"),
        legend.title = element_text(size = 6, face = "plain"),
        axis.text.x = element_text(size = 6, face = "plain"),
        plot.margin = SHARED_MARGIN)

# ----------------------------------------------------------------------------
# NAS activity column — reuse nas_activity_cascade.R (read-only) rather than
# re-deriving the same three track types on the NAS 0..8 axis. It is sourced in
# an isolated env with CASCADE_WRITE_OUTPUT=false so it builds its p1/p2/p3a/p3b
# without writing its own standalone supp file.
# ----------------------------------------------------------------------------
load_cascade <- function(script_name) {
  old_write <- Sys.getenv("CASCADE_WRITE_OUTPUT", unset = NA_character_)
  on.exit({
    if (is.na(old_write)) Sys.unsetenv("CASCADE_WRITE_OUTPUT")
    else Sys.setenv(CASCADE_WRITE_OUTPUT = old_write)
  }, add = TRUE)
  Sys.setenv(CASCADE_WRITE_OUTPUT = "false")
  env <- new.env(parent = globalenv())
  sys.source(file.path(BASE, "scripts/figures", script_name), envir = env)
  env
}
nas <- load_cascade("nas_activity_cascade.R")

# ----------------------------------------------------------------------------
# Shared legend + row-label harmonization. Both columns carry identical guide
# definitions (same labels, order, colors) so guides="collect" merges each
# legend exactly once. Row labels (horizontal, right-justified y titles) live on
# the fibrosis (left) column; the NAS (right) column drops them.
# ----------------------------------------------------------------------------
direction_labels <- c("Upregulated", "Downregulated")
nmf_breaks <- c("Inflammatory-EMT", "Fibrotic-ECM", "Hepatocyte-Metabolic",
                "Unresolved", "Noncoding", "Skeletal-muscle")
nmf_colors <- c(
  "Inflammatory-EMT" = cat_palette[3], "Fibrotic-ECM" = cat_palette[6],
  "Hepatocyte-Metabolic" = "#BDBDBD", "Unresolved" = "#D6D6D6",
  "Noncoding" = "#9E9E9E", "Skeletal-muscle" = "#7D7D7D")
celltype_breaks <- c("T cells", "Macrophages", "Fibroblasts", "Endothelial cells")
celltype_colors <- c(
  "T cells" = ct_palette[["T cells"]], "Macrophages" = ct_palette[["Macrophages"]],
  "Fibroblasts" = ct_palette[["Fibroblasts"]],
  "Endothelial cells" = ct_palette[["Endothelial cells"]])

left_y_theme <- theme(
  axis.title.y = element_text(size = 6, face = "plain", angle = 0,
                              hjust = 1, vjust = 0.5, margin = margin(r = 4)),
  plot.margin = margin(1, 2, 1, 1))
right_y_theme <- theme(axis.title.y = element_blank(),
                       plot.margin = margin(1, 2, 1, 1))
legend_theme <- theme(
  legend.position = "right", legend.justification = "top",
  legend.box = "vertical", legend.box.spacing = unit(0.5, "mm"),
  legend.spacing.y = unit(0.5, "mm"), legend.margin = margin(0, 0, 0, 1),
  legend.title = element_text(size = 6, face = "plain"),
  legend.text = element_text(size = 6, face = "plain"),
  legend.key.height = unit(2.2, "mm"), legend.key.width = unit(2.2, "mm"))

# Percent labels on the hepatocyte lines add clutter in the compact two-column
# form; drop the GeomText layers (the value axis still reads the decline).
drop_text_layers <- function(plot) {
  plot$layers <- Filter(function(l) !inherits(l$geom, "GeomText"), plot$layers)
  plot
}

# --- LEFT column: fibrosis (carries every collected legend + row labels) ------
fib_p1 <- p1 +
  scale_fill_manual(values = DEG_COLORS, breaks = levels(deg_long$direction),
                    labels = direction_labels, name = NULL) +
  labs(x = NULL, y = "DEGs") + left_y_theme
fib_p2 <- p2 +
  scale_fill_manual(values = nmf_colors, breaks = nmf_breaks,
                    name = "NMF program") +
  labs(x = NULL, y = "NMF %") + left_y_theme
fib_p3a <- drop_text_layers(p3a) +
  labs(x = NULL, y = "Hepatocytes %") + left_y_theme
fib_p3b <- p3b +
  scale_color_manual(values = celltype_colors, breaks = celltype_breaks,
                     name = "Cell type") +
  labs(x = "Fibrosis stage (F0-F4)", y = "Non-hepatocytes %") + left_y_theme

# --- RIGHT column: NAS (identical guides -> collected; no row labels) ----------
nas_p1 <- nas$p1 +
  scale_fill_manual(values = nas$deg_colors,
                    breaks = levels(nas$deg_long$direction),
                    labels = direction_labels, name = NULL) +
  labs(x = NULL, y = NULL) + guides(fill = "none") + right_y_theme
nas_p2 <- nas$p2 +
  scale_fill_manual(values = nmf_colors, breaks = nmf_breaks,
                    name = "NMF program") +
  labs(x = NULL, y = NULL) + right_y_theme
nas_p3a <- drop_text_layers(nas$p3a) +
  labs(x = NULL, y = NULL) + right_y_theme
nas_p3b <- nas$p3b +
  scale_color_manual(values = celltype_colors, breaks = celltype_breaks,
                     name = "Cell type") +
  scale_x_discrete(drop = FALSE, limits = nas$GROUPS,
                   labels = sub("^NAS", "", nas$GROUPS), expand = nas$X_EXPAND) +
  labs(x = "NAS activity group", y = NULL) + guides(color = "none") + right_y_theme

# Row-major placement keeps matched fibrosis/NAS tracks side by side.
cascade <- fib_p1 + nas_p1 +
  fib_p2 + nas_p2 +
  fib_p3a + nas_p3a +
  fib_p3b + nas_p3b +
  plot_layout(ncol = 2, widths = c(1, 1),
              heights = c(1.0, 1.0, 0.68, 0.82), guides = "collect") &
  legend_theme

out_name <- sprintf("fig3e_progression_cascade%s.pdf", file_tag)
write_output <- tolower(Sys.getenv("CASCADE_WRITE_OUTPUT", "true")) %in%
  c("true", "1", "yes")
if (write_output) {
  save_fig(cascade, file.path(PANEL_DIR, out_name),
           width = PANEL_WIDTH_IN, height = PANEL_HEIGHT_IN)

  # Underlying summary tables (caption transparency) — fibrosis column only;
  # the NAS column writes its own sidecars when nas_activity_cascade.R is run.
  fwrite(deg_counts, file.path(DATA_DIR,
    sprintf("progression_cascade_deg_counts%s.csv", file_tag)))
  fwrite(nmf_comp[order(stage, prog_label)],
         file.path(DATA_DIR, "progression_cascade_nmf_composition.csv"))
  fwrite(ct_mean, file.path(DATA_DIR, "progression_cascade_celltype_means.csv"))

  cat("Wrote", out_name, "\n")
  cat("Fibrosis DEG counts per stage (up/down):\n"); print(deg_counts)
  cat("NAS DEG counts per group (up/down):\n"); print(nas$deg_counts)
  cat(sprintf("Hepatocyte fraction F0=%.3f  F4=%.3f  (%.1f%% decline)\n",
              ct_mean[stage == "F0", Hepatocytes],
              ct_mean[stage == "F4", Hepatocytes],
              100 * (1 - ct_mean[stage == "F4", Hepatocytes] /
                         ct_mean[stage == "F0", Hepatocytes])))
}
