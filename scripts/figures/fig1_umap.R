#!/usr/bin/env Rscript
# fig1_umap.R — UMAP of integrated MASLD cohorts colored by key metadata.
#
# Uses the pre-computed Harmony-corrected UMAP from
# results/integration/umap_coordinates.csv (built by 09_batch_correction_umap.R).
# Sample IDs are aligned by row order against merged_dge.rds (1,444 raw;
# PRJNA512027 is dropped from cohort presentation, leaving 1,259 in the UMAP).
# Joins with unified_metadata.csv for fibrosis_stage / nas_score / diagnosis.
#
# Output: figures/supplementary/figS01_qc_validation/panels/fig1_umap.pdf
# (moved from main fig1; UMAP is now a supplementary QC panel.)
# Composite 2x2: (1) cohort, (2) disease state, (3) fibrosis stage, (4) diagnosis.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIGS01_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# Load + align UMAP coords with sample IDs and metadata
# ----------------------------------------------------------------------------
message("Loading integrated UMAP coordinates...")
dge <- load_merged_dge()
sample_ids <- rownames(dge$samples)
umap_dt <- fread(file.path(INT_RESULTS, "umap_coordinates.csv"))
stopifnot(nrow(umap_dt) == length(sample_ids))
umap_dt[, sample_id := sample_ids]

meta <- fread(file.path(INT_META, "unified_metadata.csv"))
plot_dt <- merge(umap_dt, meta[, .(sample_id, fibrosis_stage, nas_score,
                                   diagnosis_harmonized, age, sex)],
                 by = "sample_id", all.x = TRUE)
# umap_dt already had a `sex` column from upstream that may be missing/empty;
# prefer the metadata table's sex but coalesce.
if ("sex.x" %in% names(plot_dt) || "sex.y" %in% names(plot_dt)) {
  plot_dt[, sex := fifelse(!is.na(sex.y) & sex.y != "", sex.y, sex.x)]
  plot_dt[, c("sex.x", "sex.y") := NULL]
}

message(sprintf("Plot data (raw): %s samples across %s cohorts",
                comma(nrow(plot_dt)), uniqueN(plot_dt$dataset)))

# Drop PRJNA512027 from cohort presentation (L0/S0 library-prep batch
# perfectly confounded with diagnosis: all 34 controls L0, all 102 NASH S0).
# Still loaded by the pipeline for the fibrosis-vs-healthy contrast.
plot_dt <- plot_dt[dataset != "PRJNA512027"]
message(sprintf("Plot data (PRJNA excluded): %s samples across %s cohorts",
                comma(nrow(plot_dt)), uniqueN(plot_dt$dataset)))

# ----------------------------------------------------------------------------
# Cohort labels (short names matching Fig 1 convention)
# ----------------------------------------------------------------------------
cohort_short <- c(
  GSE126848   = "Suppli",
  GSE130970   = "Hoang",
  GSE135251   = "Govaere",
  GSE162694   = "Bril",
  GSE167523   = "Kozumi",
  GSE174478   = "Kawamura",
  GSE193066   = "Hoshida",
  GSE213621   = "Chen",
  GSE240729   = "Verschuren"
)
plot_dt[, cohort := factor(cohort_short[dataset],
                           levels = unname(cohort_short))]

# Disease state — already in `group` column from existing UMAP CSV
plot_dt[, disease_state := factor(group, levels = c("Control", "Disease"))]

# Fibrosis stage as ordered factor (F0-F4)
plot_dt[, fib_factor := factor(paste0("F", fibrosis_stage),
                               levels = paste0("F", 0:4))]
plot_dt[fibrosis_stage |> is.na(), fib_factor := NA]

# Diagnosis harmonized (Control / NAFL / Borderline / NASH)
plot_dt[, dx := factor(diagnosis_harmonized,
                       levels = c("Control", "NAFL", "Borderline", "NASH"))]

# Sex — clean blanks → NA
plot_dt[sex == "" | is.na(sex), sex := NA_character_]
plot_dt[, sex_factor := factor(sex, levels = c("F", "M"))]

# NAS score — keep numeric for gradient
plot_dt[, nas_num := suppressWarnings(as.numeric(nas_score))]

# Unified condition: use harmonized diagnosis when available; otherwise fall
# back to a fibrosis-stage tier (for cohorts annotated with fibrosis only).
# This replaces the older per-cohort lump that mixed "NAFL/NASH/Fibrosis_F*"
# into a single panel with inconsistent label semantics.
plot_dt[, unified_dx := fcase(
  diagnosis_harmonized == "Control",          "Control",
  diagnosis_harmonized == "NAFL",             "NAFL",
  diagnosis_harmonized == "Borderline",       "Borderline",
  diagnosis_harmonized == "NASH",             "NASH",
  fibrosis_stage %in% c(0, 1),                "Fibrosis F0-F1 (no NAS)",
  fibrosis_stage == 2,                        "Fibrosis F2 (no NAS)",
  fibrosis_stage %in% c(3, 4),                "Fibrosis F3-F4 (no NAS)",
  default = NA_character_
)]
plot_dt[, unified_dx := factor(unified_dx, levels = c(
  "Control", "NAFL", "Borderline", "NASH",
  "Fibrosis F0-F1 (no NAS)", "Fibrosis F2 (no NAS)", "Fibrosis F3-F4 (no NAS)"
))]

# ----------------------------------------------------------------------------
# Color palettes
# ----------------------------------------------------------------------------
# 10-cohort palette — rotate Tab10-style hues
cohort_colors <- c(
  "Suppli"     = "#1F77B4",
  "Hoang"      = "#FF7F0E",
  "Govaere"    = "#2CA02C",
  "Bril"       = "#D62728",
  "Kozumi"     = "#9467BD",
  "Kawamura"   = "#8C564B",
  "Hoshida"    = "#E377C2",
  "Chen"       = "#7F7F7F",
  "Verschuren" = "#BCBD22",
  "Gerhard"    = "#17BECF"
)

# Control is always rendered in neutral gray across all panels (project
# convention: control = gray-ish so disease/severity hues stand out).
CONTROL_GRAY <- "#9E9E9E"

disease_colors <- c(
  Control = CONTROL_GRAY,
  Disease = masld_colors$nash
)

# Fibrosis stage gradient already defined in publication_theme.R as
# fibrosis_stage_colors — F0 is the lightest (closest to control); we
# additionally tint F0 toward gray when used in this UMAP context.

dx_colors <- c(
  "Control"    = CONTROL_GRAY,
  "NAFL"       = masld_colors$nafl,
  "Borderline" = "#EC407A",
  "NASH"       = masld_colors$nash
)

sex_colors <- c(F = masld_colors$female, M = masld_colors$male)

# Unified-condition palette: Control = gray; harmonized diagnosis on a
# pink/magenta ramp; fibrosis-only fallback on a teal ramp to visually flag
# that those samples are scored on the alternate axis.
unified_dx_colors <- c(
  "Control"                  = CONTROL_GRAY,
  "NAFL"                     = "#F4A674",
  "Borderline"               = "#EC407A",
  "NASH"                     = "#C9265E",
  "Fibrosis F0-F1 (no NAS)"  = "#80CBC4",
  "Fibrosis F2 (no NAS)"     = "#26A69A",
  "Fibrosis F3-F4 (no NAS)"  = "#00695C"
)

# ----------------------------------------------------------------------------
# Helper to keep panels consistent
# ----------------------------------------------------------------------------
umap_xy <- function(p, ...) {
  p +
    coord_fixed() +
    labs(x = "UMAP 1", y = "UMAP 2", ...) +
    theme_masld(base_size = 7) +
    theme(plot.title = element_text(size = 8, face = "bold"),
          legend.title = element_text(size = 6.5, face = "bold"),
          legend.text  = element_text(size = 6),
          legend.key.size = unit(0.28, "cm"))
}

dot <- function(aes_call, color_override = NULL, ...) {
  layer <- if (is.null(color_override))
    geom_point(aes_call, size = 0.4, alpha = 0.7, shape = 16)
  else
    geom_point(aes_call, size = 0.4, alpha = 0.7, shape = 16,
               color = color_override)
  rasterize_layer(layer)
}

# ----------------------------------------------------------------------------
# Panel 1 — by cohort
# ----------------------------------------------------------------------------
p1 <- umap_xy(
  ggplot(plot_dt, aes(UMAP1, UMAP2, color = cohort)) +
    rasterize_layer(geom_point(size = 0.35, alpha = 0.7, shape = 16)) +
    scale_color_manual(values = cohort_colors, name = "Cohort", drop = FALSE) +
    guides(color = guide_legend(ncol = 2,
                                override.aes = list(size = 1.4, alpha = 1))),
  title = "Cohort")

# ----------------------------------------------------------------------------
# Panel 2 — by disease state
# ----------------------------------------------------------------------------
p2 <- umap_xy(
  ggplot(plot_dt, aes(UMAP1, UMAP2, color = disease_state)) +
    rasterize_layer(geom_point(size = 0.35, alpha = 0.7, shape = 16)) +
    scale_color_manual(values = disease_colors, name = "State",
                       na.value = "gray85") +
    guides(color = guide_legend(override.aes = list(size = 1.4, alpha = 1))),
  title = "Disease state")

# ----------------------------------------------------------------------------
# Panel 3 — by fibrosis stage (1,149 samples; rest gray)
# ----------------------------------------------------------------------------
fib_dt <- plot_dt[order(!is.na(fib_factor))]   # NAs first → drawn beneath
p3 <- umap_xy(
  ggplot(fib_dt, aes(UMAP1, UMAP2, color = fib_factor)) +
    rasterize_layer(geom_point(size = 0.35, alpha = 0.75, shape = 16)) +
    scale_color_manual(values = fibrosis_stage_colors, name = "Fibrosis",
                       na.value = "gray85", drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 1.4, alpha = 1))),
  title = sprintf("Fibrosis stage (n = %s)",
                  comma(sum(!is.na(plot_dt$fib_factor)))))

# ----------------------------------------------------------------------------
# Panel 4 — by diagnosis (Control / NAFL / Borderline / NASH)
# ----------------------------------------------------------------------------
dx_dt <- plot_dt[order(!is.na(dx))]
p4 <- umap_xy(
  ggplot(dx_dt, aes(UMAP1, UMAP2, color = dx)) +
    rasterize_layer(geom_point(size = 0.35, alpha = 0.75, shape = 16)) +
    scale_color_manual(values = dx_colors, name = "Diagnosis",
                       na.value = "gray85", drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 1.4, alpha = 1))),
  title = sprintf("Diagnosis (n = %s)",
                  comma(sum(!is.na(plot_dt$dx)))))

# ----------------------------------------------------------------------------
# Panel 5 — by sex
# ----------------------------------------------------------------------------
sex_dt <- plot_dt[order(!is.na(sex_factor))]
p5 <- umap_xy(
  ggplot(sex_dt, aes(UMAP1, UMAP2, color = sex_factor)) +
    rasterize_layer(geom_point(size = 0.35, alpha = 0.75, shape = 16)) +
    scale_color_manual(values = sex_colors, name = "Sex",
                       na.value = "gray85", drop = FALSE,
                       labels = c(F = "Female", M = "Male")) +
    guides(color = guide_legend(override.aes = list(size = 1.4, alpha = 1))),
  title = sprintf("Sex (n = %s)",
                  comma(sum(!is.na(plot_dt$sex_factor)))))

# ----------------------------------------------------------------------------
# Panel 6 — by NAS score (continuous gradient)
# ----------------------------------------------------------------------------
nas_dt <- plot_dt[order(!is.na(nas_num))]
p6 <- umap_xy(
  ggplot(nas_dt, aes(UMAP1, UMAP2, color = nas_num)) +
    rasterize_layer(geom_point(size = 0.35, alpha = 0.75, shape = 16)) +
    scale_color_gradientn(
      colors = c("#FFF8E1", "#FFD54F", "#F57F17", "#C9265E", "#A01753"),
      limits = c(0, 8), breaks = c(0, 2, 4, 6, 8),
      name = "NAS",
      na.value = "gray85"),
  title = sprintf("NAS score (n = %s)",
                  comma(sum(!is.na(plot_dt$nas_num)))))

# ----------------------------------------------------------------------------
# Panel 7 — unified condition (harmonized diagnosis + fibrosis-only fallback)
# Replaces the older condition_lump panel that mixed NAFL/NASH with raw
# Fibrosis_F* labels inconsistently. See unified_dx classifier above.
# ----------------------------------------------------------------------------
cond_dt <- plot_dt[order(!is.na(unified_dx))]
p7 <- umap_xy(
  ggplot(cond_dt, aes(UMAP1, UMAP2, color = unified_dx)) +
    rasterize_layer(geom_point(size = 0.35, alpha = 0.75, shape = 16)) +
    scale_color_manual(values = unified_dx_colors,
                       name = "Condition", na.value = "gray85",
                       drop = FALSE) +
    guides(color = guide_legend(ncol = 1,
                                override.aes = list(size = 1.4, alpha = 1))),
  title = sprintf("Unified condition (n = %s)",
                  comma(sum(!is.na(plot_dt$unified_dx)))))

# ----------------------------------------------------------------------------
# Compose 2 x 4 (8th cell empty)
# ----------------------------------------------------------------------------
spacer <- patchwork::plot_spacer()
fig <- (p1 | p7 | p2 | p5) / (p4 | p3 | p6 | spacer) +
  plot_annotation(
    title    = "Integrated atlas — UMAP of 1,259 QC-passing samples",
    subtitle = "Harmony batch-corrected on dataset; 9 cohorts (PRJNA512027 excluded)",
    theme = theme(plot.title    = element_text(size = 9, face = "bold",
                                               family = "Helvetica"),
                  plot.subtitle = element_text(size = 7, color = "gray35",
                                               family = "Helvetica"))
  )

save_fig(fig, file.path(PANEL_DIR, "fig1_umap.pdf"),
         width = fig_full_width * 1.6, height = 4.6)

fp <- file.path(PANEL_DIR, "fig1_umap.pdf")
if (file.exists(fp)) {
  message(sprintf("\nOutput: %s (%s)", fp,
                  utils:::format.object_size(file.size(fp), "auto")))
}

# Also dump per-panel CSV so the metadata join is reproducible
fwrite(plot_dt[, .(sample_id, dataset, cohort, UMAP1, UMAP2,
                   disease_state, condition, unified_dx,
                   fibrosis_stage, fib_factor,
                   nas_score, diagnosis_harmonized, dx, age,
                   sex, sex_factor)],
       file.path(PANEL_DIR, "fig1_umap_data.csv"))
message("Wrote fig1_umap_data.csv")
