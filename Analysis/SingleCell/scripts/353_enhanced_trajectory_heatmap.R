#!/usr/bin/env Rscript
# ============================================================================
# 353_enhanced_trajectory_heatmap.R
#
# Upgraded main stage-trajectory CCC heatmap with proper validation overlays.
# Replaces the cramped 3-track right strip in 349 / 349b with a 7-track
# annotation block, built around the dedup top-30 INDEPENDENT signals.
#
# Panel A : main heatmap (rows = top-30 dedup independent signals,
#           cols = 4 disease stages, fill = donor-mean z-scored
#           -log10 magnitude_rank, RdBu diverging midpoint=0 capped at +/-2)
# Panel A-side : 7 validation tracks, each in its own column block
#           (>=0.4 inch wide), readable not cramped:
#             1. Coarse SH beta (RdBu)
#             2. Pseudotime beta (RdBu)
#             3. F-stage beta (RdBu)
#             4. Bulk concordance (4 stage axes, tick marks)
#             5. LOO replication rate (viridis 0-1)
#             6. Bootstrap CI excludes zero (binary)
#             7. Dissoc status (categorical green/yellow/red/gray)
# Panel B : per-row stacked summary "concordant / tested axes"
# Panel C : legend block
#
# Output:
#   figures/supplementary/stage_ccc/figS_stage_ccc_trajectory_v2.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(viridis)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")
SUPP_DIR <- file.path(BASE, "figures/supplementary/stage_ccc")
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)

theme_src <- file.path(BASE, "scripts/figures/publication_theme.R")
if (file.exists(theme_src)) source(theme_src)

cap <- function(x, lim) pmin(pmax(x, -lim), lim)

# MASLD project palette (publication_theme.R `masld_colors`)
DIV_LOW   <- masld_colors$down       # "#1565C0" blue
DIV_HIGH  <- masld_colors$up         # "#C2185B" magenta
OK_GREEN  <- masld_colors$conserved  # "#00695C" teal (validates)
WARN_MAG  <- masld_colors$up         # magenta (caveat)
NEUTRAL   <- masld_colors$ns         # "#BDBDBD" gray

# ---------------------------------------------------------------------------
# Load all inputs
# ---------------------------------------------------------------------------
DEDUP    <- fread(file.path(OUT_DIR, "dedup_top30_independent_signals.tsv"))
LR_LONG  <- fread(file.path(OUT_DIR, "all_donor_lr_scores.tsv.gz"))
META_EXT <- fread(file.path(OUT_DIR, "donor_metadata_extended.tsv"))
COARSE   <- fread(file.path(OUT_DIR, "stage_lr_lmm_coarse.tsv"))
CONTIN   <- fread(file.path(OUT_DIR, "stage_lr_lmm_continuous.tsv"))
FSTAGE   <- fread(file.path(OUT_DIR, "stage_lr_lmm_fstage.tsv"))
BOOT     <- fread(file.path(OUT_DIR, "bootstrap_ci.tsv"))
DISSOC   <- fread(file.path(OUT_DIR, "all_dissoc_sensitivity.tsv"))
LOO      <- fread(file.path(OUT_DIR, "loo_replication_rate_per_lr.tsv"))
BULK_CC  <- fread(file.path(OUT_DIR, "lr_bulk_concordance.tsv"))

cat(sprintf("[353] loaded %d dedup rows\n", nrow(DEDUP)))

# Build canonical row id
DEDUP[, lr_pair := paste(representative_ligand, receptor, sep = "__")]
DEDUP[, fam_tag := ifelse(n_members >= 2,
                          sprintf(" (+%d %s)", n_members - 1, ligand_family),
                          "")]
DEDUP[, row_label := sprintf("%s | %s%s",
                             ct_pair,
                             paste(representative_ligand, receptor, sep = "->"),
                             fam_tag)]

# ---------------------------------------------------------------------------
# Heatmap matrix: donor-mean z-scored -log10(magnitude_rank) per 4 stages
# ---------------------------------------------------------------------------
LR_LONG[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
LR_LONG[, ct_pair := paste(source, target, sep = "->")]
LR_LONG <- merge(LR_LONG,
                 META_EXT[, .(sample, disease_stage_coarse)],
                 by = "sample", all.x = TRUE)
LR_LONG[, score := -log10(pmax(magnitude_rank, 1e-4))]

top_keys <- unique(DEDUP[, .(ct_pair, lr_pair)])
hm <- merge(LR_LONG, top_keys, by = c("ct_pair", "lr_pair"))
hm_mean <- hm[, .(score_mean = mean(score, na.rm = TRUE),
                  n_donors   = .N),
              by = .(ct_pair, lr_pair, disease_stage_coarse)]
hm_mean[, z := scale(score_mean)[, 1], by = .(ct_pair, lr_pair)]
hm_mean[, disease_stage_coarse := factor(disease_stage_coarse,
        levels = c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"))]
hm_mean <- merge(hm_mean,
                 DEDUP[, .(ct_pair, lr_pair, row_label)],
                 by = c("ct_pair", "lr_pair"))

# Row ordering: sort by SH coarse beta (descending)
DEDUP_ORD <- DEDUP[order(-Estimate)]
row_levels <- DEDUP_ORD$row_label

# ---------------------------------------------------------------------------
# Build annotation table per row
# ---------------------------------------------------------------------------
ann <- DEDUP[, .(ct_pair, lr_pair, row_label)]
ann[, coarse_SH := DEDUP$Estimate]

# Pseudotime continuous beta
CONTIN_SH <- CONTIN[term == "macrophage_pseudotime_mean",
                    .(ct_pair = paste(source, target, sep = "->"),
                      lr_pair = paste(ligand_complex, receptor_complex, sep = "__"),
                      cont_beta = Estimate)]
ann <- merge(ann, CONTIN_SH, by = c("ct_pair", "lr_pair"), all.x = TRUE)

# F-stage beta
FS <- FSTAGE[term == "F_stage_numeric",
             .(ct_pair = paste(source, target, sep = "->"),
               lr_pair = paste(ligand_complex, receptor_complex, sep = "__"),
               fstage_beta = Estimate)]
ann <- merge(ann, FS, by = c("ct_pair", "lr_pair"), all.x = TRUE)

# Bulk concordance per axis -- 4 axes
bulk_long <- unique(
  BULK_CC[, .(ct_pair = paste(source, target, sep = "->"),
              lr_pair = paste(ligand_complex, receptor_complex, sep = "__"),
              axis, both_concordant)]
)
ann_keys <- unique(ann[, .(ct_pair, lr_pair)])
bulk_long <- merge(bulk_long, ann_keys, by = c("ct_pair", "lr_pair"))
# Map axis -> short label
axis_short <- c(
  "Steatosis_vs_Healthy"            = "Steat",
  "Steatohepatitis_vs_Healthy"      = "SH",
  "Cirrhosis_vs_Healthy"            = "Cirr",
  "macrophage_pseudotime_all_MASLD" = "AllMASLD"
)
bulk_long[, axis_short := axis_short[axis]]
bulk_long <- bulk_long[!is.na(axis_short)]
# For each (row, axis), TRUE if ANY entry concordant (the table is per direction)
bulk_wide <- dcast(bulk_long,
                   ct_pair + lr_pair ~ axis_short,
                   value.var = "both_concordant",
                   fun.aggregate = function(x) any(as.logical(x), na.rm = TRUE),
                   fill = FALSE)
# Ensure all 4 cols
for (s in c("Steat", "SH", "Cirr", "AllMASLD")) {
  if (!s %in% names(bulk_wide)) bulk_wide[[s]] <- FALSE
}
ann <- merge(ann, bulk_wide, by = c("ct_pair", "lr_pair"), all.x = TRUE)
for (s in c("Steat", "SH", "Cirr", "AllMASLD")) {
  ann[[s]][is.na(ann[[s]])] <- FALSE
}

# LOO replication rate (strict)
LOO_SUB <- LOO[, .(ct_pair, lr_pair,
                   loo_rate = replication_rate_strict)]
ann <- merge(ann, LOO_SUB, by = c("ct_pair", "lr_pair"), all.x = TRUE)

# Bootstrap CI excludes zero
BOOT_SUB <- BOOT[, .(ct_pair, lr_pair, boot_ok = as.logical(ci_excludes_zero))]
ann <- merge(ann, BOOT_SUB, by = c("ct_pair", "lr_pair"), all.x = TRUE)

# Dissoc status: take the SH-vs-Healthy classification first; fallback to Cirr
DISSOC[, dissoc_class := classification]
dissoc_sh <- DISSOC[contrast == "SH_vs_Healthy",
                    .(ct_pair, lr_pair,
                      dissoc_sh = dissoc_class)]
dissoc_ci <- DISSOC[contrast == "Cirrhosis_vs_Healthy",
                    .(ct_pair, lr_pair,
                      dissoc_ci = dissoc_class)]
ann <- merge(ann, dissoc_sh, by = c("ct_pair", "lr_pair"), all.x = TRUE)
ann <- merge(ann, dissoc_ci, by = c("ct_pair", "lr_pair"), all.x = TRUE)
ann[, dissoc := fifelse(!is.na(dissoc_sh), dissoc_sh,
                  fifelse(!is.na(dissoc_ci), dissoc_ci, NA_character_))]
ann[, dissoc := fifelse(is.na(dissoc), "not_tested", dissoc)]

# ---------------------------------------------------------------------------
# Warning flag & star
# ---------------------------------------------------------------------------
ann[, boot_ok_filled := fifelse(is.na(boot_ok), FALSE, boot_ok)]
ann[, warn := (dissoc %in% c("weakened", "reversed")) |
              (!boot_ok_filled)]
ann[, star := !is.na(loo_rate) & loo_rate >= 1.0 - 1e-6]

# Compute concordance summary across testable axes
# Tested axes: pseudotime (have value), fstage (have value),
# bulk SH (TRUE/FALSE; FALSE counts as tested-but-failed only if entry exists -- we treat all 4 bulk as tested),
# LOO (always tested), bootstrap (always tested if in BOOT), dissoc (tested if not "not_tested")
concord_summary <- function(r) {
  tested <- 0L; concord <- 0L
  # SH coarse beta -- always tested, concordant if matches expected direction
  # We treat SH coarse beta as the anchor (always counted: concordant if |Estimate|>0)
  # but to avoid double-counting against itself, skip.
  if (!is.na(r$cont_beta)) {
    tested <- tested + 1L
    if (sign(r$cont_beta) == sign(r$coarse_SH)) concord <- concord + 1L
  }
  if (!is.na(r$fstage_beta)) {
    tested <- tested + 1L
    if (sign(r$fstage_beta) == sign(r$coarse_SH)) concord <- concord + 1L
  }
  # Bulk axes: each axis is one test
  for (s in c("Steat", "SH", "Cirr", "AllMASLD")) {
    if (s %in% names(r)) {
      tested <- tested + 1L
      if (isTRUE(r[[s]])) concord <- concord + 1L
    }
  }
  if (!is.na(r$loo_rate)) {
    tested <- tested + 1L
    if (r$loo_rate >= 5/6 - 1e-6) concord <- concord + 1L
  }
  if (!is.na(r$boot_ok)) {
    tested <- tested + 1L
    if (isTRUE(r$boot_ok)) concord <- concord + 1L
  }
  if (r$dissoc != "not_tested") {
    tested <- tested + 1L
    if (r$dissoc == "survives_dissoc") concord <- concord + 1L
  }
  list(concord = concord, tested = tested)
}

cs <- lapply(seq_len(nrow(ann)), function(i) concord_summary(as.list(ann[i])))
ann[, concord_n := sapply(cs, function(x) x$concord)]
ann[, tested_n  := sapply(cs, function(x) x$tested)]

# Convergent definition (per spec):
# significant on >=2 axes (coarse SH + pseudotime or F-stage)
# AND bootstrap survives AND dissoc-robust AND LOO >= 5/6
# Helper: NA-safe logical
isTRUE_vec <- function(x) !is.na(x) & as.logical(x)

ann[, sig_axes :=
       as.integer(!is.na(coarse_SH)) +
       as.integer(!is.na(cont_beta) & sign(cont_beta) == sign(coarse_SH)) +
       as.integer(!is.na(fstage_beta) & sign(fstage_beta) == sign(coarse_SH))]
ann[, convergent :=
       sig_axes >= 2 &
       isTRUE_vec(boot_ok) &
       dissoc == "survives_dissoc" &
       !is.na(loo_rate) & loo_rate >= 5/6 - 1e-6]

# ---------------------------------------------------------------------------
# Final row labels with star prefix + red color flag
# ---------------------------------------------------------------------------
ann[, row_label_disp := ifelse(star,
                               paste0("★ ", row_label),
                               row_label)]
# Map original row_label -> display label
disp_lookup <- setNames(ann$row_label_disp, ann$row_label)
row_levels_disp <- disp_lookup[row_levels]

hm_mean[, row_label := factor(disp_lookup[row_label], levels = row_levels_disp)]
ann[, row_label := factor(row_label_disp, levels = row_levels_disp)]

# Row text colour
ann[, label_color := ifelse(warn, WARN_MAG, "grey20")]
row_colors <- setNames(ann$label_color[match(row_levels_disp, ann$row_label)],
                       row_levels_disp)

# ---------------------------------------------------------------------------
# PANEL A: main heatmap
# ---------------------------------------------------------------------------
p_hm <- ggplot(hm_mean,
               aes(x = disease_stage_coarse, y = row_label, fill = cap(z, 2))) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HIGH,
                       midpoint = 0,
                       limits = c(-2, 2), oob = scales::squish,
                       name = "z(-log10 mag.rank)") +
  scale_x_discrete(labels = c("Healthy", "Steat.", "SH", "Cirr.")) +
  labs(x = NULL, y = NULL,
       title = sprintf("Stage-trajectory CCC: top %d independent signals",
                       nrow(DEDUP))) +
  theme_masld(base_size = 7) +
  theme(axis.text.y      = element_text(size = 5.5,
                                        color = row_colors[row_levels_disp]),
        axis.text.x      = element_text(angle = 35, hjust = 1, size = 6.5),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        legend.position  = "top",
        legend.key.width = unit(0.65, "cm"),
        legend.key.height= unit(0.22, "cm"),
        legend.title     = element_text(size = 6, face = "bold"),
        legend.text      = element_text(size = 5.5),
        legend.margin    = margin(0, 0, 0, 0),
        plot.title       = element_text(face = "bold", size = 8, hjust = 0,
                                        margin = margin(b = 2)),
        plot.title.position = "plot",
        plot.margin      = margin(2, 2, 2, 2))

# ---------------------------------------------------------------------------
# Generic single-track helper for diverging RdBu betas
# ---------------------------------------------------------------------------
make_beta_track <- function(df, col, title, cap_at = 2,
                            show_y = FALSE, legend_pos = "none") {
  d <- df[, c("row_label", col), with = FALSE]
  setnames(d, col, "val")
  d[, val := cap(val, cap_at)]
  ggplot(d, aes(x = 1, y = row_label, fill = val)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HIGH,
                         midpoint = 0,
                         limits = c(-cap_at, cap_at),
                         oob = scales::squish,
                         na.value = "grey90", name = "β") +
    scale_x_continuous(breaks = 1, labels = title, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_masld(base_size = 6) +
    theme(axis.text.y     = if (show_y)
                              element_text(size = 6,
                                           color = row_colors[row_levels_disp])
                            else element_blank(),
          axis.text.x     = element_text(angle = 35, hjust = 1, size = 5.5),
          axis.ticks      = element_blank(),
          axis.line       = element_blank(),
          legend.position = legend_pos,
          legend.key.width  = unit(0.4, "cm"),
          legend.key.height = unit(0.2, "cm"),
          legend.title    = element_text(size = 5.5, face = "bold"),
          legend.text     = element_text(size = 5),
          plot.margin     = margin(2, 1, 2, 1))
}

p_t1 <- make_beta_track(ann, "coarse_SH",  "SH beta",       cap_at = 2)
p_t2 <- make_beta_track(ann, "cont_beta",  "Pseudotime",    cap_at = 2)
p_t3 <- make_beta_track(ann, "fstage_beta","F-stage",       cap_at = 1)

# Bulk concordance: 4 stage cols as binary ticks
bulk_m <- melt(ann[, .(row_label, Steat, SH, Cirr, AllMASLD)],
               id.vars = "row_label",
               variable.name = "bulk_axis", value.name = "ok")
bulk_m[, ok := as.logical(ok)]
bulk_m[, bulk_axis := factor(bulk_axis,
                             levels = c("Steat", "SH", "Cirr", "AllMASLD"))]

p_t4 <- ggplot(bulk_m, aes(x = bulk_axis, y = row_label, fill = ok)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_manual(values = c("TRUE"  = OK_GREEN,
                               "FALSE" = "grey95"),
                    na.value = "grey90",
                    name = "Bulk concord.") +
  geom_text(data = bulk_m[ok == TRUE],
            aes(label = "✓"), size = 1.9, color = "white",
            fontface = "bold") +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.y     = element_blank(),
        axis.text.x     = element_text(angle = 35, hjust = 1, size = 5.5),
        axis.ticks      = element_blank(),
        axis.line       = element_blank(),
        legend.position = "none",
        plot.margin     = margin(2, 1, 2, 1))

# LOO replication rate (mako sequential — same family as panel C in 351)
p_t5 <- ggplot(ann, aes(x = 1, y = row_label, fill = loo_rate)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1),
                       na.value = "grey90", name = "LOO rep.") +
  geom_text(aes(label = ifelse(is.na(loo_rate), "",
                               sprintf("%.2f", loo_rate))),
            size = 1.7, color = "white") +
  scale_x_continuous(breaks = 1, labels = "LOO rep.", expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.y     = element_blank(),
        axis.text.x     = element_text(angle = 35, hjust = 1, size = 5.5),
        axis.ticks      = element_blank(),
        axis.line       = element_blank(),
        legend.position = "none",
        plot.margin     = margin(2, 1, 2, 1))

# Bootstrap CI excludes zero
p_t6 <- ggplot(ann, aes(x = 1, y = row_label, fill = boot_ok)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_manual(values = c("TRUE"  = masld_colors$control,  # deep blue
                               "FALSE" = "grey95"),
                    na.value = "grey90",
                    name = "Boot CI") +
  geom_text(data = ann[!is.na(boot_ok) & boot_ok == TRUE],
            aes(label = "✓"), size = 1.9, color = "white",
            fontface = "bold") +
  scale_x_continuous(breaks = 1, labels = "Boot CI", expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.y     = element_blank(),
        axis.text.x     = element_text(angle = 35, hjust = 1, size = 5.5),
        axis.ticks      = element_blank(),
        axis.line       = element_blank(),
        legend.position = "none",
        plot.margin     = margin(2, 1, 2, 1))

# Dissoc status categorical (green = survives; magenta gradient = caveats)
dissoc_colors <- c("survives_dissoc" = OK_GREEN,         # teal: passes
                   "weakened"        = masld_colors$masl, # pale magenta
                   "reversed"        = masld_colors$up,   # magenta
                   "not_tested"      = "grey85")
p_t7 <- ggplot(ann, aes(x = 1, y = row_label, fill = dissoc)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_manual(values = dissoc_colors, name = "Dissoc.") +
  scale_x_continuous(breaks = 1, labels = "Dissoc.", expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.y     = element_blank(),
        axis.text.x     = element_text(angle = 35, hjust = 1, size = 5.5),
        axis.ticks      = element_blank(),
        axis.line       = element_blank(),
        legend.position = "none",
        plot.margin     = margin(2, 1, 2, 1))

# ---------------------------------------------------------------------------
# PANEL B: per-row stacked summary "concordant / tested"
# ---------------------------------------------------------------------------
bar_df <- ann[, .(row_label,
                  concord = concord_n,
                  remain  = tested_n - concord_n)]
bar_long <- melt(bar_df, id.vars = "row_label",
                 variable.name = "kind", value.name = "count")
bar_long[, kind := factor(kind, levels = c("concord", "remain"),
                          labels = c("concordant", "tested but not"))]

p_bar <- ggplot(bar_long, aes(x = count, y = row_label, fill = kind)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.1) +
  geom_text(data = ann,
            aes(x = tested_n + 0.3, y = row_label,
                label = sprintf("%d/%d", concord_n, tested_n)),
            inherit.aes = FALSE,
            size = 1.9, hjust = 0, color = "grey25") +
  scale_fill_manual(values = c("concordant"     = OK_GREEN,
                               "tested but not" = "grey80")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)),
                     breaks = pretty_breaks(3)) +
  labs(x = "axes concordant", y = NULL, fill = NULL,
       title = "Validation summary") +
  theme_masld(base_size = 6) +
  theme(axis.text.y        = element_blank(),
        axis.text.x        = element_text(size = 5.5),
        axis.title.x       = element_text(size = 6),
        axis.ticks.y       = element_blank(),
        axis.line.y        = element_blank(),
        legend.position    = "top",
        legend.key.size    = unit(0.25, "cm"),
        legend.text        = element_text(size = 5.5),
        legend.margin      = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.1, "cm"),
        plot.title         = element_text(face = "bold", size = 7, hjust = 0,
                                          margin = margin(b = 2)),
        plot.title.position = "plot",
        plot.margin        = margin(2, 2, 2, 1))

# ---------------------------------------------------------------------------
# PANEL C: legend block (a free-form ggplot with no data; text annotation)
# ---------------------------------------------------------------------------
legend_lines <- c(
  "Heatmap fill: z-scored -log10(magnitude_rank) per LR pair across donors,",
  "  capped at ±2. Red = stronger interaction; blue = weaker.",
  "",
  "Annotation tracks (left to right):",
  "  1) Coarse SH β    -- lmer beta, Steatohepatitis vs Healthy",
  "  2) Pseudotime β   -- lmer beta on macrophage pseudotime (continuous)",
  "  3) F-stage β     -- lmer beta on documented F-stage 0-4",
  "  4) Bulk concord.  -- per-stage RNA-seq LFC sign-match (4 axes)",
  "  5) LOO rep.       -- leave-one-dataset-out replication rate (0-1, strict)",
  "  6) Boot CI        -- 95% bootstrap CI on SH-β excludes zero",
  "  7) Dissoc.        -- dissociation-bias sensitivity classification",
  "",
  "Symbols:",
  "  ★ row prefix      = LOO replication 6/6 (perfect)",
  "  red row label    = dissoc flagged (weakened/reversed) OR boot CI crosses 0",
  "  green / yellow / red = survives / weakened / reversed under dissoc",
  "  gray             = not tested",
  "",
  "Convergent (per panel B): ≥2 sig axes + boot survives + dissoc-robust + LOO≥5/6"
)
p_legend <- ggplot() +
  annotate("text", x = 0, y = seq_along(legend_lines),
           label = rev(legend_lines),
           hjust = 0, size = 1.8, family = "sans", color = "grey25") +
  scale_y_continuous(limits = c(0, length(legend_lines) + 1)) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(2, 4, 2, 4))

# ---------------------------------------------------------------------------
# Assemble layout
#
# Top row:  [ p_hm ][ p_t1 ][ p_t2 ][ p_t3 ][ p_t4 ][ p_t5 ][ p_t6 ][ p_t7 ][ p_bar ]
# Bottom row: legend block spans full width
# ---------------------------------------------------------------------------
top_row <- p_hm +
  p_t1 + p_t2 + p_t3 + p_t4 + p_t5 + p_t6 + p_t7 + p_bar +
  plot_layout(nrow = 1,
              widths = c(3.6, 0.3, 0.3, 0.3, 1.2, 0.4, 0.3, 0.45, 1.4))

full <- top_row / p_legend +
  plot_layout(heights = c(14, 2)) +
  plot_annotation(
    title    = "Stage-trajectory CCC: top 30 signals with multi-axis validation",
    subtitle = "Rows ordered by SH-vs-Healthy coarse β (descending)",
    theme    = theme(plot.title    = element_text(face = "bold", size = 9,
                                                  hjust = 0,
                                                  margin = margin(b = 1)),
                     plot.subtitle = element_text(size = 6.5, color = "grey30",
                                                  hjust = 0,
                                                  margin = margin(b = 2)),
                     plot.title.position = "plot")
  )

out_pdf <- file.path(SUPP_DIR, "figS_stage_ccc_trajectory_v2.pdf")
ggsave(out_pdf, full, width = 10, height = 8, device = cairo_pdf)
cat(sprintf("[output] %s\n", out_pdf))

# ---------------------------------------------------------------------------
# Side TSV: convergence audit for the report
# ---------------------------------------------------------------------------
audit <- ann[, .(row_label, ct_pair, lr_pair,
                 coarse_SH, cont_beta, fstage_beta,
                 Steat, SH, Cirr, AllMASLD,
                 loo_rate, boot_ok, dissoc,
                 sig_axes, concord_n, tested_n, convergent, star)]
audit_tsv <- file.path(OUT_DIR, "top30_validation_audit.tsv")
fwrite(audit, audit_tsv, sep = "\t")
cat(sprintf("[output] %s\n", audit_tsv))

# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
cat("\n========== TOP-30 VALIDATION SUMMARY ==========\n")
cat(sprintf("Total rows                        : %d\n", nrow(ann)))
cat(sprintf("Convergent (strict, all 4 tests)  : %d\n", sum(ann$convergent)))
cat(sprintf("Marginal (failing >=1 test)       : %d\n", sum(!ann$convergent)))
cat("\nBreakdown:\n")
cat(sprintf("  bootstrap CI excludes zero      : %d / %d\n",
            sum(ann$boot_ok, na.rm = TRUE),
            sum(!is.na(ann$boot_ok))))
cat(sprintf("  dissoc survives                 : %d / %d (tested)\n",
            sum(ann$dissoc == "survives_dissoc"),
            sum(ann$dissoc != "not_tested")))
cat(sprintf("  dissoc weakened                 : %d\n",
            sum(ann$dissoc == "weakened")))
cat(sprintf("  dissoc reversed                 : %d\n",
            sum(ann$dissoc == "reversed")))
cat(sprintf("  LOO >=5/6 strict                : %d / %d\n",
            sum(!is.na(ann$loo_rate) & ann$loo_rate >= 5/6 - 1e-6),
            sum(!is.na(ann$loo_rate))))
cat(sprintf("  LOO 6/6 perfect (star)          : %d\n", sum(ann$star)))
cat(sprintf("  Sig axes >=2 (coarse+contin/F)  : %d\n",
            sum(ann$sig_axes >= 2)))
cat("===============================================\n")
cat("\n[done]\n")
